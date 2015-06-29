package com.butlitsky.mk.ensembles;

import com.butlitsky.mk.IEnsemble;
import com.butlitsky.mk.options.CLOptions;
import org.apache.commons.math3.util.FastMath;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.text.NumberFormat;
import java.util.List;

import static org.apache.commons.math3.util.FastMath.ceil;

/**
 * Writing and loading new particles configurations and other particles related info
 * (correlations etc.).
 * <p/>
 * Also responsible for initial coordinates setting.
 * <p/>
 * Created by aristofun on 24.09.14.
 */
public class GibbsConfigurationManager {
    private NumberFormat short_format, long_format;

    private final String myFolder, configFile;
    private Path myConfigPath;
    // private Path myCorrPath;


    private BufferedWriter longTailWriter;
    // writer for current point 2 box parameters file
    private BufferedWriter additionalStateWriter;

    private final int Nei;
    private final double[][][] prtcls;
    private boolean saveLongTail = false;

    private final GibbsEnsemble myEnsemble;

// -–––– correlation stuff --------------- // todo: fix correlation later
//    private final double[] corrNormirovka = new double[2]; // 0 - first box, 1 - second box
//    private final double[] corrDr = new double[2];
//    private final double[][][] corrArray = new double[2][3][];
//    private int corrAverager = 0;


    public GibbsConfigurationManager(GibbsEnsemble ensemble) {
        myEnsemble = ensemble;

        myFolder = ensemble.myFolder;
        configFile = IEnsemble.STATE_FILE;

        Nei = ensemble.opt.getNumParticles() / 2;

        myConfigPath = GibbsConfigurationManager.getPath(myFolder + "/" + configFile);

        short_format = myEnsemble.SHORT_FORMAT;
        long_format = myEnsemble.FORMAT;

        prtcls = myEnsemble.getParticles();

        // OPTIONS first bit == save longtail
        saveLongTail = ((myEnsemble.opt.getStrategy() & 1) == 1);

//        corrNormirovka = 4. * (boxSize * boxSize * boxSize) / (Nei * Nei);
//        corrDr = Math.sqrt(3.0 * boxSize * boxSize) / 1.99999999999 / (CORR_LENGTH);
//        corrArray[0] = new double[CORR_LENGTH];
//        corrArray[1] = new double[CORR_LENGTH];
//        corrArray[2] = new double[CORR_LENGTH];

    }

    /**
     * reads config or generates new random if any errors ocurred
     *
     * @throws java.io.IOException on errors reading file
     */
    void loadConfiguration() throws IOException {
//        myCorrPath = EnsembleController.getPath(myFolder + "/" + CORR_FILE);
        Files.createDirectories(GibbsConfigurationManager.getPath(myFolder));

        boolean recover = myEnsemble.opt.isOld();

        // reading old state
        if (recover) {
            if (Files.exists(myConfigPath)) {
                try {
                    readCoordinates(Files.readAllLines(myConfigPath, Charset.forName("UTF-8")));
                    System.out.print(": config loaded - ok");
                } catch (Exception e) {
                    System.out.println("WARNING: failed to read config for " + myFolder);
                    System.out.println(e.getLocalizedMessage());
                    recover = false;
                }
            } else {
                recover = false;
            }
        }

        if (!recover) {
            System.out.print("! from scratch");
            initParticlesPosition(); // initial distribution, reset counters
            saveConfiguration();
//            System.exit(0);
        }
    }

    /**
     * Loads state from given strings array, treating first line as a general parameters
     * and all other as an array of particles' coordinates inside correspondent boxes.
     * <p/>
     * throws Exception if no valid config could be read (either
     */
    private void readCoordinates(List<String> strings) throws Exception {
//          first line format:
// current step, boxBorder, avg. energy 1 (per prtcl), avg density 1, avg. energy 2, avg. density 2, total Gamma
        String firstline = strings.remove(0);

        myEnsemble.setCurrStep(Integer.parseInt(firstline.split("\\s+")[0]));
        final int boxBorder = Integer.parseInt(firstline.split("\\s+")[1]);
        myEnsemble.setBoxBorder(boxBorder);

        final double avgE1 = Double.parseDouble(firstline.split("\\s+")[2]);
        final double density1 = Double.parseDouble(firstline.split("\\s+")[3]);
        final double avgE2 = Double.parseDouble(firstline.split("\\s+")[4]);
        final double density2 = Double.parseDouble(firstline.split("\\s+")[5]);

        myEnsemble.setCurrReducedEnergies(avgE1, avgE2);
        myEnsemble.setDensitiesAvg(density1, density2);

        if (strings.size() != Nei * 2) {
            throw new IndexOutOfBoundsException("file size doesn't fit particles number");
        }

        for (int type = 0; type < 2; type++) {
            for (int i = 0 + type * Nei; i < Nei + type * Nei; i++) {
                String[] strs = strings.get(i).split("\\s");

                prtcls[type][0][i - type * Nei] = Double.parseDouble(strs[0]);
                prtcls[type][1][i - type * Nei] = Double.parseDouble(strs[1]);
                prtcls[type][2][i - type * Nei] = Double.parseDouble(strs[2]);
            }
        }
    }

    /**
     * fill random arrays of particles, resetting all counters
     */
    void initParticlesPosition() {
        final double[] boxSizes = myEnsemble.getBoxSizes();
        final int boxBorder = myEnsemble.getBoxBorder();

        // uniformly distribute particles in the box NaCL structure (fcc)
        if (CLOptions.START_FROM_FCC) {
            fillBoxFCC(0, boxBorder, boxSizes[0]);
            fillBoxFCC(boxBorder, Nei, boxSizes[1]);
        } else {
            // random configuration.
            for (int type = 0; type < 2; type++) {
                for (int j = 0; j < Nei; j++) {
                    int box = (j < boxBorder) ? 0 : 1;
                    /* spreading the particles */
                    prtcls[type][0][j] = myEnsemble.myRandom(boxSizes[box]);
                    prtcls[type][1][j] = myEnsemble.myRandom(boxSizes[box]);
                    prtcls[type][2][j] = myEnsemble.myRandom(boxSizes[box]);
                }
            }
        }
    }

    private void fillBoxFCC(int start, int stop, double box_size) {
        final int length = stop - start;
        final int Nx = (int) ceil(FastMath.cbrt(length * 2));
        final double Lx = box_size / Nx;

        for (int j = 0; j < length * 2; j++) {
            int zLayer = j / (Nx * Nx);
            int yLayer = (j % (Nx * Nx)) / Nx;
            int xLayer = (j % (Nx * Nx)) % Nx;

            int type = ((j % 2) == 0) ? 0 : 1;  // ion, electron, ion, electron...
            prtcls[type][0][start + j / 2] = xLayer * Lx;
            prtcls[type][1][start + j / 2] = yLayer * Lx;
            prtcls[type][2][start + j / 2] = zLayer * Lx;

        }
    }

    void saveConfiguration() {
        // create strings list from arrays
        try {
            BufferedWriter writer = Files.newBufferedWriter(myConfigPath, Charset.forName("UTF-8"));
            final double avgEnergy1 = myEnsemble.getAvgEnergy(0);
            final double avgEnergy2 = myEnsemble.getAvgEnergy(1);
//          first line format:
//
// current step, boxBorder, avg. energy 1 (per prtcl), avg density 1, avg. energy 2, avg. density 2, total Gamma
            writer.write(
                    "" + myEnsemble.getCurrStep() + "\t"
                            + myEnsemble.getBoxBorder() + "\t"
                            + long_format.format(avgEnergy1) + "\t"
                            + long_format.format(myEnsemble.getDensitiesAvg()[0]) + "\t"
                            + long_format.format(avgEnergy2) + "\t"
                            + long_format.format(myEnsemble.getDensitiesAvg()[1]) + "\t"
                            + short_format.format(myEnsemble.opt.getGamma()) + "\t"
            );

            writer.newLine();
            writeCoordinates(writer);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: failed to save state for " + myFolder);
        }
    }

    /**
     * Writing current coordinates to BufferedWriter
     *
     * @param writer
     * @throws IOException
     */
    private void writeCoordinates(BufferedWriter writer) throws IOException {
        // all E goes first, then all Ions (NOT box after box!)
        for (int type = 0; type < 2; type++) {
            for (int i = 0 + type * Nei; i < Nei + type * Nei; i++) {
                writer.write(
                        long_format.format(prtcls[type][0][i - type * Nei]) + "\t"
                                + long_format.format(prtcls[type][1][i - type * Nei]) + "\t"
                                + long_format.format(prtcls[type][2][i - type * Nei])
                );
                writer.newLine();
            }
        }
    }

    void applyAdditionalStrategies() {
//        if (saveLongTail) {
/*            try {
                longTailWriter = Files.newBufferedWriter(
                        GibbsConfigurationManager.getPath(myFolder + "/" + IEnsemble.LONGTAIL_FILE),
                        Charset.defaultCharset(), StandardOpenOption.CREATE,
                        opt.isOld() ? StandardOpenOption.APPEND : StandardOpenOption.TRUNCATE_EXISTING,
                        StandardOpenOption.WRITE);
            } catch (Exception e) {
                System.out.println("ERROR: Can't create " + LONGTAIL_FILE + " for " + myFolder);
            }*/
//        }

        // open current point adiitonal parameters logging "2box.dat"
        try {
            additionalStateWriter = Files.newBufferedWriter(
                    GibbsConfigurationManager.getPath(myFolder + "/" + GibbsEnsemble.GIBBS_STATE_FILE),
                    Charset.forName("UTF-8"), StandardOpenOption.CREATE,
                    myEnsemble.opt.isOld() ? StandardOpenOption.APPEND : StandardOpenOption.TRUNCATE_EXISTING,
                    StandardOpenOption.WRITE);
        } catch (Exception e) {
            System.out.println("ERROR: Can't create " + GibbsEnsemble.GIBBS_STATE_FILE + " for " + myFolder);
        }
    }

    void closeLongTails() {
/*        if (saveLongTail) try {
            longTailWriter.flush();
            longTailWriter.close();
        } catch (IOException e1) {
            System.out.println("ERROR: can't close long tail writer for " + myFolder);
        } */

        // closing gibbs 2 box state file
        try {
            additionalStateWriter.flush();
            additionalStateWriter.close();
        } catch (IOException e) {
            System.out.println("ERROR: can't close " + myFolder + " gibbs state writer!");
            e.printStackTrace();
        }
    }

    void workOnMidCalc() {
/*        if (saveLongTail) try {
            writeCoordinates(longTailWriter);
            longTailWriter.flush();
        } catch (IOException e1) {
            System.out.println("ERROR: can't write " + myFolder + " long tail to writer!");
        } */
    }


    public void workOnFrequentCalc(int curr_step) {
        // write reduced densities and particles numbers in the following format
        // curr_step, N1, v*1, N2, v*2

        try {
            additionalStateWriter.write(
                    curr_step + "\t" +
                            myEnsemble.getBoxBorder() * 2 + "\t" +
                            long_format.format(myEnsemble.getVstar(0)) + "\t" +
                            (Nei - myEnsemble.getBoxBorder()) * 2 + "\t" +
                            long_format.format(myEnsemble.getVstar(1)) + "\t"
            );
            additionalStateWriter.newLine();
//            System.out.println("written");
        } catch (IOException e) {
            System.out.println("ERROR: can't write " + myFolder + " gibbs state to " +
                                       GibbsEnsemble.GIBBS_STATE_FILE);
            e.printStackTrace();
        }
    }


    /**
     * Utility method for file seeking
     */
    public static final Path getPath(String path) {
        return FileSystems.getDefault().getPath(".", path);
    }

}


// TODO: fix correlation stuff

        /*private final void fillCorrAray(int i, int j, int corrIndex) {
        if (i != j) {
            // get the index of a radius in array
            final int idx =
                    (int) (Math.sqrt(
                            dSquared(
                                    Xs[i] - Xs[j], Ys[i] - Ys[j],
                                    Zs[i] - Zs[j])) / corrDr);
            // increment ION-ION array
            if (idx < CORR_LENGTH)
                corrArray[corrIndex][idx]++;
            else
                System.out.println("WARNING: " + myFolder + " correlation index out of bounds!");
        }
    } */

    /*
    private final void calcCorrelation() {
        // get the simple corrArray for ION-ION
        for (int i = 0; i < (Nei / 2); i++) {
            for (int j = 0; j < (Nei / 2); j++) {
                fillCorrAray(i, j, 0);
            }
        }
        // get the simple corrArray for ELECTRON-ELECTRON
        for (int i = (Nei / 2); i < Nei; i++) {
            for (int j = (Nei / 2); j < Nei; j++) {
                fillCorrAray(i, j, 2);
            }
        }

        // get the simple corrArray for ELECTRON-ION
        for (int i = (Nei / 2); i < Nei; i++) {
            for (int j = 0; j < (Nei / 2); j++) {
                fillCorrAray(i, j, 1);
            }
        }
        // increment the number of correlation array calculatings
        corrAverager++;
    }
    */

 /*   private final void saveCorrelation() {
        List<String> strings = new ArrayList<>(CORR_LENGTH);

        if (corrAverager == 0) {
            System.out.println(
                    "WARNING: corrAverager == 0 for " + myFolder + ", " +
                            "skip saveCorrelation");
            return;
        }

        for (int i = 0; i < CORR_LENGTH; i++) {
            // radius in the middle of a sherical layer
            final double r = i * corrDr + 0.5 * corrDr;
//                                              #    sherical layer volume     #
            final double norm = corrNormirovka / (4.0 * Math.PI * r * r * corrDr * corrAverager);
            // writing
            strings.add(
                    FORMAT.format(r) + "\t"
                            + FORMAT.format(corrArray[0][i] * norm) + "\t"
                            + FORMAT.format(corrArray[1][i] * norm) + "\t"
                            + FORMAT.format(corrArray[2][i] * norm) + "\t"
            );
        }

        try {
            Files.write(myCorrPath, strings, Charset.defaultCharset());
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: failed to save correlation for " + myFolder);
        }
    }
 */