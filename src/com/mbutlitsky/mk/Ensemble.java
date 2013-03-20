package com.mbutlitsky.mk;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static com.mbutlitsky.mk.EnsembleController.getPath;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

/**
 * Created with IntelliJ IDEA.
 * User: aristofun
 * Date: 03.03.13
 * Time: 13:04
 */
public abstract class Ensemble implements IEnsemble {
    private static final int SAVE_CONFIG_INT = 501;
    private static final int SAVE_CORR_INT = 1200;
    private static final int CALC_CORR_INT = 600;
    public static double DELTA_FACTOR = 1.3;
    /**
     * overrides specific particle number in ini file if set in CLI options
     */
    public static int DEFAULT_NUM_PARTICLES = 0;

    private final NumberFormat FORMAT = new DecimalFormat(EOptions.SCIENTIFIC_FORMAT_STR);
    private final NumberFormat SHORT_FORMAT = new DecimalFormat(EOptions.SHORT_FORMAT_STR);

    // ----------  controllers  --------------
    private final EOptions opt;
    private final String myFolder;
    private final MersenneTwisterFast rnd;
    private Path myConfigPath;
    private Path myCorrPath;
    private volatile boolean finished = false;


    // ------------ Monte Karlo --------------
    private final int numPart;
    private final int numSteps;
    private int currStep;
    private int which;
    private double xTrial;
    private double yTrial;
    private double zTrial;

    // ------------ mathematics --------------
    protected final int T;
    private final double boxSize;
    private final double halfBox;
    private final double maxDelta;
    private final double avgDistance;

    private final double corrNormirovka;
    private final double corrDr;
    private final double[][] corrArray = new double[3][];
    private int corrAverager = 0;

    private final double[] Xs;
    private final double[] Ys;
    private final double[] Zs;
    private double energy;


    public Ensemble(EOptions options) {
        rnd = new MersenneTwisterFast();
        opt = options;
        boxSize = pow(opt.getNumParticles() / (2 * opt.getDensity()), 0.3333333333333333) / BOHR;
        halfBox = boxSize / 2.0;
        avgDistance = pow(opt.getDensity(), -0.333333333333333333) / BOHR;

        // ignore specific delta factor if set to zero
        maxDelta = (opt.getMaxDelta() == 0.0) ? DELTA_FACTOR : opt.getMaxDelta();

        numPart = (DEFAULT_NUM_PARTICLES == 0) ? opt.getNumParticles() : DEFAULT_NUM_PARTICLES;
        numSteps = opt.getNumSteps();
        T = opt.getT();
        myFolder = opt.getFolder();

//        double maxDistance = sqrt(3.0 * boxSize * boxSize) / 1.99;

        corrNormirovka = 4. * (boxSize * boxSize * boxSize) / (numPart * numPart);
        corrDr = sqrt(3.0 * boxSize * boxSize) / 1.99999999999 / (CORR_LENGTH);

        corrArray[0] = new double[CORR_LENGTH];
        corrArray[1] = new double[CORR_LENGTH];
        corrArray[2] = new double[CORR_LENGTH];

        // initiating particles configuration
        Xs = new double[numPart];
        Ys = new double[numPart];
        Zs = new double[numPart];

        initFromStateFile();
        resetEnergy();
    }

    public static final void fillArray(double[] dabls, String line) {
        String[] strs = line.split("\\s");
        for (int i = 0; i < dabls.length; i++) {
            dabls[i] = Double.parseDouble(strs[i]);
        }
    }

    protected final String doubleArrayToString(double[] dabls) {
        String result = "";
        for (double dabl : dabls) {
            result += FORMAT.format(dabl) + "\t";
        }
        return result;
    }

    /**
     * reads config or generates new random if any errors ocurred
     */
    private void initFromStateFile() {
        myConfigPath = getPath(myFolder + "/" + STATE_FILE);
        myCorrPath = getPath(myFolder + "/" + CORR_FILE);

        try {
            Files.createDirectories(getPath(myFolder));
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: can't find " + myFolder + " directory, " +
                    "worker is dead");
            finished = true;
            return;
        }

        // reading old state
        if (opt.isOld()) {
            if (Files.exists(myConfigPath)) {
                try {
                    initArrays(Files.readAllLines(myConfigPath, Charset.defaultCharset()));
                } catch (Exception e) {
                    System.out.println("WARNING: failed to read config for " + myFolder);
                    System.out.println(e.getLocalizedMessage());
                    opt.setOld(false);
                }
            } else {
                opt.setOld(false);
                System.out.println("WARNING: initial config not found for " + myFolder + ", " +
                        "starting from scratch");
            }
        }

        if (!opt.isOld()) initRandomState(); // initial distribution, reset counters
    }

    /**
     * fill random arrays of particles, resetting all counters
     */
    private void initRandomState() {
        currStep = 0;

        // random configuration.
        for (int j = 0; j < numPart; j++) {
            /* spreading the particles */
            Xs[j] = myRandom(boxSize);
            Ys[j] = myRandom(boxSize);
            Zs[j] = myRandom(boxSize);
        }
    }

    /**
     * calculates full system energy from scratch, based on current particles configuration
     */
    private final void resetEnergy() {
        energy = 0;
        for (int i = 0; i < numPart; i++) {
            for (int j = i + 1; j < numPart; j++) {
                energy += getPotential(i, j);
            }
        }
    }

    /**
     * trialX/Y/Z must belong to j particle
     *
     * @return potential between two particles assuming j[trialX, trialY, trialZ] particle
     */
    private final double getPotential(int i, int j, double trialX, double trialY, double trialZ) {
        final double r = sqrt(dSquared(trialX - Xs[i], trialY - Ys[i], trialZ - Zs[i]));

        if (i != j) {
            if (j < (numPart / 2))   // First _num/2 are IONS
            {
                if (i < (numPart / 2)) // ION-ION
                    return getPotential(r, false);
                else              // ION - electron
                    return getPotential(r, true);
            } else                   // Last _num/2 are Electrons
            {
                if (i < (numPart / 2)) // Electron - ION
                    return getPotential(r, true);
                else              // Electron - Electron
                    return getPotential(r, false);
            }
        }
        return 0;
    }

    /**
     * returns current energy between to particles defined by i,j indexes
     * if i == j returns 0
     */
    private final double getPotential(int i, int j) {
        return getPotential(i, j, Xs[j], Ys[j], Zs[j]);
    }

    protected abstract double getPotential(double r, boolean attraction);

    public final double dSquared(double dx, double dy, double dz) {
        dx = dSquared(dx);
        dy = dSquared(dy);
        dz = dSquared(dz);
        return ((dx * dx) + (dy * dy) + (dz * dz));
    }

    private final double dSquared(double dx) {
        dx = Math.abs(dx);
        return (dx > halfBox) ? (halfBox - dx % halfBox) : dx;
    }

    /**
     * random (-1.0; 1.0)
     */
    protected final double myRandom() {
        return ((double) rnd.nextInt() / Integer.MAX_VALUE);
//        return ThreadLocalRandom.current().nextDouble(-1.,1.);
    }

    protected final int myRandom(int size) {
        return rnd.nextInt(size);
//        return ThreadLocalRandom.current().nextInt(size);
    }

    /**
     * returns Mersenne Twister random (0; size) double value
     */
    protected final double myRandom(double size) {
        // faster on Core 2 Duo
        return size * ((double) rnd.nextInt() / Integer.MAX_VALUE + 1) / 2;
//        return rnd.nextDouble() * size;
//        return ThreadLocalRandom.current().nextDouble(size);
    }

    /**
     * throws IOException if no valid config could be read
     */
    private void initArrays(List<String> strings) throws Exception {
        // first line: curr_step energy
        String step_en = strings.remove(0);
        currStep = Integer.parseInt(step_en.split("\\s+")[0]);

        fillArray(Xs, strings.get(0));
        fillArray(Ys, strings.get(1));
        fillArray(Zs, strings.get(2));
    }

    private final void saveCorrelation() {
        List<String> strings = new ArrayList<String>(CORR_LENGTH);

        if (corrAverager == 0) {
            System.out.println("WARNING: corrAverager == 0 for " + myFolder + ", " +
                    "skip saveCorrelation");
            return;
        }

        for (int i = 0; i < CORR_LENGTH; i++) {
            // radius in the middle of a sherical layer
            final double r = i * corrDr + 0.5 * corrDr;
//                                              #    sherical layer volume     #
            final double norm = corrNormirovka / (4.0 * Math.PI * r * r * corrDr * corrAverager);
            // writing
            strings.add(FORMAT.format(r) + "\t"
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

    public void saveState() {
        resetEnergy();
        // create strings list from arrays
        try {
            Files.write(myConfigPath, flushArrays(), Charset.defaultCharset());
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: failed to save state for " + myFolder);
        }
    }

    private final List<String> flushArrays() {
        return Arrays.asList(
                "" + currStep + "\t" + SHORT_FORMAT.format(energy / numPart) + "\t"
                        + SHORT_FORMAT.format(opt.getGamma()),
                doubleArrayToString(Xs),
                doubleArrayToString(Ys),
                doubleArrayToString(Zs)
        );
    }

    @Override
    public final int getCurrStep() {
        return currStep;  //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * new blank thread is running this method until stopped or comes to the numSteps
     * It quietly stops on any exceptions with message to system out
     */
    @Override
    public void run() {
        for (int i = currStep; i < numSteps; i++) {
            if (finished) {
                saveState();
                System.out.println("Exitting by stop " + myFolder + ", 'finished' found 'true'");
                return;
            }

            if (i % CALC_CORR_INT == 0) calcCorrelation();
            if (i % SAVE_CORR_INT == 0) saveCorrelation();
            if (i % SAVE_CONFIG_INT == 0) {
                saveState();
            }

            play();
            currStep = i;
        }

        finished = true;
        saveState();
        System.out.println("Ensemble " + myFolder + " finished after " + numSteps + " steps.");
    }

    private final void fillCorrAray(int i, int j, int corrIndex) {
        if (i != j) {
            // get the index of a radius in array
            final int idx =
                    (int) (sqrt(dSquared(Xs[i] - Xs[j], Ys[i] - Ys[j], Zs[i] - Zs[j])) / corrDr);
            // increment ION-ION array
            if (idx < CORR_LENGTH)
                corrArray[corrIndex][idx]++;
            else
                System.out.println("WARNING: " + myFolder + " correlation index out of bounds!");
        }
    }

    private final void calcCorrelation() {
        // get the simple corrArray for ION-ION
        for (int i = 0; i < (numPart / 2); i++) {
            for (int j = 0; j < (numPart / 2); j++) {
                fillCorrAray(i, j, 0);
            }
        }
        // get the simple corrArray for ELECTRON-ELECTRON
        for (int i = (numPart / 2); i < numPart; i++) {
            for (int j = (numPart / 2); j < numPart; j++) {
                fillCorrAray(i, j, 2);
            }
        }

        // get the simple corrArray for ELECTRON-ION
        for (int i = (numPart / 2); i < numPart; i++) {
            for (int j = 0; j < (numPart / 2); j++) {
                fillCorrAray(i, j, 1);
            }
        }
        // increment the number of correlation array calculatings
        corrAverager++;
    }

    /**
     * performs one act of monte-karlo play randomly moves each particle and saves the state
     * with some probability
     */
    private final void play() {
        for (int j = 0; j < numPart; j++) {
            // move random particle
            final double deltaE = moveParticle();
            // transition probability checking
            // energy increased
            if (deltaE > 0) {
                // compare the transition probability with random
                // All energies are in kT
                if (Math.exp(-deltaE) >= myRandom(1.0))   // accept movement
                    acceptTrial(deltaE);

                // energy decreased, accept configuration
            } else {
                acceptTrial(deltaE);
            }
        }
    }

    private final void acceptTrial(double deltaE) {
        Xs[which] = xTrial;
        Ys[which] = yTrial;
        Zs[which] = zTrial;
        energy += deltaE;
    }

    /**
     * returns trial energy shift for moved particle
     */
    private final double moveParticle() {
        final double x = myRandom() * maxDelta;
        final double y = myRandom() * maxDelta;
        final double z = myRandom() * maxDelta;

        // setting new trial index and coordinates
        which = myRandom(numPart);
        xTrial = correctPosition(Xs[which] + x);
        yTrial = correctPosition(Ys[which] + y);
        zTrial = correctPosition(Zs[which] + z);

        // Calculating the energy shift
        double oldE = 0, newE = 0;

        // firstly, old Energy
        for (int i = 0; i < numPart; i++) {
            oldE += getPotential(i, which);
        }

        // then, new Energy
        for (int i = 0; i < numPart; i++) {
            newE += getPotential(i, which, xTrial, yTrial, zTrial);
        }
        return newE - oldE;
    }

    /**
     * set right position inside the Box
     */
    private final double correctPosition(double coord) {
        return ((boxSize + coord % boxSize) % boxSize);
        // optimize it
        // int mod (int a, int b)
//        {
//            int ret = a % b;
//            if(ret < 0)
//                ret+=b;
//            return ret;
//        }
// OLD code:   (coord - boxSize * (long) (coord * 2 / boxSize - 1));
    }

    @Override
    public double getEnergy() {
        return energy;
    }

    @Override
    public double getPressure() {
        return 0;
    }

    @Override
    public final boolean isFinished() {
        return finished;
    }

    @Override
    public final EOptions getOptions() {
        return opt;
    }

    @Override
    public final void stop() {
        finished = true;
    }

    public String toString() {
        return "ensmbl:" + myFolder;
    }
}
