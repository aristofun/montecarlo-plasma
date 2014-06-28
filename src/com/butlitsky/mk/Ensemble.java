package com.butlitsky.mk;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;
import java.util.concurrent.ThreadLocalRandom;

import static com.butlitsky.mk.EnsembleController.getPath;
import static org.apache.commons.math3.util.FastMath.exp;
import static org.apache.commons.math3.util.FastMath.pow;

/**
 * total num steps usually
 * 5 mln (500 particles x 10000 steps) to 10 mln (1000 particles x 10000 stps).
 * <p/>
 * User: aristofun
 * Date: 03.03.13
 * Time: 13:04
 */
public abstract class Ensemble implements IEnsemble {

    /**
     * scaling factor to convert energy value to kT units
     * SCALE_FACTOR == e^2 / (Bohr * k)
     */
    public static final double SCALE_FACTOR = 315775.01611746440408;
    //    private static final int CALC_ENERGY_INT = 72073; // make less for big num part
    private final int SAVE_CONFIG_N_CORR_INT;// = 12287; // make less for big num part
    //    private static final int CALC_CORR_n_ENERGY_INT = 91;
    private final int CALC_CORR_INT;// = 361275;
    //    private static final int CALC_ENERGY_INT = 12251;
    private final int AVERAGE_ENERGY_INT;// = 1831;

    private static final int INITIAL_NUM_STEPS = 30;  // how many steps from beginning to ignore


    protected final NumberFormat FORMAT = new DecimalFormat(EOptions.SCIENTIFIC_FORMAT_STR);
    protected final NumberFormat SHORT_FORMAT = new DecimalFormat(EOptions.SHORT_FORMAT_STR);

    // ----------  controllers  --------------
    private final EOptions opt;
    private final String myFolder;

    private final MersenneTwisterFast rnd;
    private final ThreadLocalRandom localRnd;

    private Path myConfigPath;
    private Path myCorrPath;
    private volatile boolean finished = false;
    private final boolean saveLongTail;
    private BufferedWriter longTailWriter;


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

    private final double corrNormirovka;
    private final double corrDr;
    private final double[][] corrArray = new double[3][];
    private int corrAverager = 0;

    protected final double[] Xs;
    protected final double[] Ys;
    protected final double[] Zs;

//    private double potential = 0;

    public static int NUM_ENERGY_AVG_POINTS = -1;

    // averaging energies values stack
    private final Deque<Double> energies;
    private double avgEnergy = 0;
    private double energy = 0;


    public Ensemble(EOptions options) {
        rnd = new MersenneTwisterFast();
        localRnd = ThreadLocalRandom.current();

        opt = options;

        numPart = opt.getNumParticles();
        numSteps = opt.getNumSteps() * numPart;
        T = opt.getT();

        boxSize = pow(numPart / (2 * opt.getDensity()),
                0.33333333333333333333) / BOHR; // parameter is Ne(Ni),
        // we double 'cause total density is twice bigger
        halfBox = boxSize / 2.0;

        // average e-i distance => x 2
        double avgDistance = pow(2 * opt.getDensity(), -0.3333333333333333333333) / BOHR;

        // ignore specific delta factor if set to zero

        double factor = opt.getMaxDelta();

        maxDelta = (factor == 0.0) ? boxSize : ((factor >= 1) ? factor * avgDistance : factor * boxSize);

        myFolder = opt.getFolder();

        System.out.print(myFolder + ": avg dist = " + SHORT_FORMAT.format(avgDistance) + ", " +
                "delta = " + SHORT_FORMAT.format(maxDelta) + ", boxSize = "
                + SHORT_FORMAT.format(boxSize));

        corrNormirovka = 4. * (boxSize * boxSize * boxSize) / (numPart * numPart);
        corrDr = StrictMath.sqrt(3.0 * boxSize * boxSize) / 1.99999999999 / (CORR_LENGTH);

        corrArray[0] = new double[CORR_LENGTH];
        corrArray[1] = new double[CORR_LENGTH];
        corrArray[2] = new double[CORR_LENGTH];

        // initiating particles configuration
        Xs = new double[numPart];
        Ys = new double[numPart];
        Zs = new double[numPart];

        // OPTIONS first bit == save longtail
        saveLongTail = ((options.getStrategy() & 1) == 1);
        // simulation timeframes scaled by number of particles
        CALC_CORR_INT = numPart;
        AVERAGE_ENERGY_INT = numPart * 3;
        SAVE_CONFIG_N_CORR_INT = numPart * 5;

        if (NUM_ENERGY_AVG_POINTS < 0) {
            NUM_ENERGY_AVG_POINTS = numSteps - INITIAL_NUM_STEPS * numPart;
        }
        System.out.print(", AVG_POINTS: " + NUM_ENERGY_AVG_POINTS);
        energies = new ArrayDeque<>(NUM_ENERGY_AVG_POINTS);
    }

    /**
     * Must be run by children after super constructor. Default strategy – set up initial state from file
     * or random if failed to read config from disk.
     */
    protected void initialize() {
        loadFromStateFile();
        initEnergy();

        applyAdditionalStrategies();
    }

    /**
     * sets average potential from file or resets it to zero
     */
    private void initEnergy() {

        if (avgEnergy == 0) {
            calcEnergy();
            averageEnergy();
        } else {
            energy = avgEnergy;
        }
    }


    private void applyAdditionalStrategies() {
        if (saveLongTail) {
            try {
                longTailWriter = Files.newBufferedWriter(getPath(myFolder + "/" + LONGTAIL_FILE),
                        Charset.defaultCharset(), StandardOpenOption.CREATE,
                        opt.isOld() ? StandardOpenOption.APPEND : StandardOpenOption.TRUNCATE_EXISTING,
                        StandardOpenOption.WRITE);
            } catch (Exception e) {
                System.out.println("ERROR: Can't create " + LONGTAIL_FILE + " for " + myFolder);
            }
        }
    }

    /**
     * reads config or generates new random if any errors ocurred
     */
    private void loadFromStateFile() {
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
                    loadArrays(Files.readAllLines(myConfigPath, Charset.defaultCharset()));
                } catch (Exception e) {
                    System.out.println("WARNING: failed to read config for " + myFolder);
                    System.out.println(e.getLocalizedMessage());
                    opt.setOld(false);
                }
            } else {
                opt.setOld(false);
                System.out.print("! from scratch: " + myFolder + "\t ");
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
     * Record current energy value to averager array
     */
    private final void calcEnergy() {
        energy = getCurrentEnergy();
        oldEnergy();
    }

    private final void oldEnergy() {
        if (energies.size() > NUM_ENERGY_AVG_POINTS - 1) energies.pollLast();
        energies.addFirst(energy);
    }

    /**
     * calculates full system potential from scratch, based on current particles configuration
     */
    private final void averageEnergy() {
        double enrg = 0.0;

        for (Double val : energies) {
            enrg = enrg + val;
        }

        avgEnergy = enrg / energies.size();
    }

    protected double getCurrentEnergy() {
        double newEn = 0;
        for (int i = 0; i < numPart; i++) {
            for (int j = i + 1; j < numPart; j++) {
                newEn = newEn + getEnergy(i, j);
            }
        }
        return newEn;
    }

    /**
     * trialX/Y/Z must belong to j particle
     *
     * @return potential between two particles assuming j[trialX, trialY, trialZ] particle
     */
    private final double getPotential(int i, int j, double trialX, double trialY, double trialZ) {
        final double r = StrictMath.sqrt(dSquared(trialX - Xs[i], trialY - Ys[i], trialZ - Zs[i]));

        if (i != j) {
            if (j < (numPart / 2))   // First _num/2 are IONS
            {
                if (i < (numPart / 2)) // ION-ION
                    return getPotentialAsym(r, false, true);
                else              // ION - electron
                    return getPotentialAsym(r, false, false);
            } else                   // Last _num/2 are Electrons
            {
                if (i < (numPart / 2)) // Electron - ION
                    return getPotentialAsym(r, false, false);
                else              // Electron - Electron
                    return getPotentialAsym(r, true, false);
            }
        }
        return 0;
    }

    /**
     * returns current potential between to particles defined by i,j indexes
     * if i == j returns 0
     */
    private final double getPotential(int i, int j) {
        return getPotential(i, j, Xs[j], Ys[j], Zs[j]);
    }

    protected abstract double getPotential(double r, boolean attraction);

    protected abstract double getPotentialAsym(double r, boolean ee, boolean ii);

    private final double getEnergy(int i, int j) {
        final double r = StrictMath.sqrt(dSquared(Xs[j] - Xs[i], Ys[j] - Ys[i], Zs[j] - Zs[i]));

        if (i != j) {
            if (j < (numPart / 2))   // First _num/2 are IONS
            {
                if (i < (numPart / 2)) // ION-ION
                    return getEnergyAsym(r, false, true);
                else              // ION - electron
                    return getEnergyAsym(r, false, false);
            } else                   // Last _num/2 are Electrons
            {
                if (i < (numPart / 2)) // Electron - ION
                    return getEnergyAsym(r, false, false);
                else              // Electron - Electron
                    return getEnergyAsym(r, true, false);
            }
        }
        return 0;
    }

    protected abstract double getEnergy(double r, boolean attraction);

    protected abstract double getEnergyAsym(double r, boolean ee, boolean ii);

    protected final double dSquared(double dx, double dy, double dz) {
        dx = fit2box(dx);
        dy = fit2box(dy);
        dz = fit2box(dz);
        return ((dx * dx) + (dy * dy) + (dz * dz));
    }

    protected final double fit2box(double dx) {
//        dx = abs(dx);
//        return (dx > halfBox) ? (halfBox - dx % halfBox) : dx;
        if (dx > halfBox) {
            dx -= boxSize;
        } else if (dx < -halfBox) {
            dx += boxSize;
        }

        return dx;
    }

    /**
     * random (-1.0; 1.0)
     */
    final double myRandom() {
        return rnd.nextDouble(true, true) * 2.0 - 1.0;
//        return ((double) rnd.nextInt() / Integer.MAX_VALUE);
//        return ThreadLocalRandom.current().nextDouble(-1.,1.);
    }

    /**
     * returns Mersenne Twister random [0; size) double value
     */
    final double myRandom(double size) {
        return rnd.nextDouble() * size;

//      faster on Core 2 Duo
//        return size * ((double) rnd.nextInt() / Integer.MAX_VALUE + 1) / 2;
    }

    /**
     * throws IOException if no valid config could be read
     */
    private void loadArrays(List<String> strings) throws Exception {
        // first line: curr_step potential
        String step_en = strings.remove(0);
        // #step, average potential, ...
        currStep = Integer.parseInt(step_en.split("\\s+")[0]);
        avgEnergy = Double.parseDouble(step_en.split("\\s+")[1]);

        if (strings.size() != numPart) {
            throw new Exception("file size doesn't fit particles number");
        }

        for (int i = 0; i < strings.size(); i++) {
            String[] strs = strings.get(i).split("\\s");
            Xs[i] = Double.parseDouble(strs[0]);    //
            Ys[i] = Double.parseDouble(strs[1]);    //
            Zs[i] = Double.parseDouble(strs[2]);    //
        }
    }

    private final void saveCorrelation() {
        List<String> strings = new ArrayList<>(CORR_LENGTH);

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

    void saveState() {
        averageEnergy();

        // create strings list from arrays
        try {
            BufferedWriter writer = Files.newBufferedWriter(myConfigPath, Charset.defaultCharset());

            writer.write("" + currStep + "\t"
                    + FORMAT.format(avgEnergy) + "\t"
                    + SHORT_FORMAT.format(avgEnergy / numPart) + "\t"
                    + SHORT_FORMAT.format(opt.getGamma()));

            writer.newLine();
            writeStateTo(writer);
            writer.close();
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: failed to save state for " + myFolder);
        }
    }

    private void writeStateTo(BufferedWriter writer) throws IOException {
        for (int i = 0; i < numPart; i++) {
            writer.write(FORMAT.format(Xs[i]) + "\t"
                    + FORMAT.format(Ys[i]) + "\t" + FORMAT.format(Zs[i]));
            writer.newLine();
        }
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
        System.out.println();
        if (currStep >= numSteps || finished) {
            System.out.print(myFolder + " No run.\t");
            finished = true;
            return;
        }

        int i = currStep;

        // fast run through initial steps, unstable configuration
        if (currStep < INITIAL_NUM_STEPS * numPart && numSteps > INITIAL_NUM_STEPS * numPart) {
            System.out.println("Ignoring first " + INITIAL_NUM_STEPS * numPart + " steps...");

            while (i < INITIAL_NUM_STEPS * numPart) {
                if (finished) {
                    System.out.print("STOP " + myFolder + ", finished=true\t");
                    break;
                }
                play();
                i++;
            }
        }

        currStep = i;

        if (!finished) {
            while (i < numSteps) {
                if (finished) {
                    System.out.print("STOP " + myFolder + ", finished=true\t");
                    break;
                }

                if (play()) {
                    calcEnergy();
                    currStep = i;
                } else {
                    oldEnergy();
                }

                if (i % CALC_CORR_INT == 0) { // no need for correlation accuracy
                    calcCorrelation();
                }


                if (i % AVERAGE_ENERGY_INT == 0) {
                    averageEnergy();
                    if (saveLongTail) saveLongTail();
                }

                if (i % SAVE_CONFIG_N_CORR_INT == 0) {
                    saveCorrelation();
                    saveState();
                }

                i++;
            }

        }
        currStep = i;
        finished = true;

        saveState();
        saveLongTail();
        closeLongTail();

        System.out.print("" + myFolder + " finished on " + currStep + " steps.\t");
    }

    private void closeLongTail() {
        if (saveLongTail) try {
            longTailWriter.flush();
            longTailWriter.close();
        } catch (IOException e1) {
            System.out.println("ERROR: can't close long tail writer for " + myFolder);
        }
    }

    private void saveLongTail() {
        try {
            if (saveLongTail) {
                writeStateTo(longTailWriter);
                longTailWriter.flush();
            }
        } catch (IOException e1) {
            System.out.println("ERROR: can't write " + myFolder + " long tail to writer!");
        }
    }

    private final void fillCorrAray(int i, int j, int corrIndex) {
        if (i != j) {
            // get the index of a radius in array
            final int idx =
                    (int) (StrictMath.sqrt(dSquared(Xs[i] - Xs[j], Ys[i] - Ys[j],
                            Zs[i] - Zs[j])) / corrDr);
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
     * performs one act of monte-karlo play randomly moves one particle and saves the state
     * with some probability
     */
    private final boolean play() {
//        for (int j = 0; j < numPart; j++) {
        // move random particle
        final double deltaE = moveParticle();
        // transition probability checking
        // potential increased
        if (deltaE > 0) {
            // compare the transition probability with random
            // All energies are in kT
            if (exp(-deltaE) >= localRnd.nextDouble(1.0)) {
                acceptTrial();
                return true;
            }
            // potential decreased, accept configuration
        } else {
            acceptTrial();
            return true;
        }

        return false;
    }

    private final void acceptTrial() {
        Xs[which] = xTrial;
        Ys[which] = yTrial;
        Zs[which] = zTrial;
    }

    /**
     * returns trial potential shift for moved particle
     */
    private final double moveParticle() {
        final double x = myRandom() * maxDelta;
        final double y = myRandom() * maxDelta;
        final double z = myRandom() * maxDelta;

        // setting new trial index and coordinates
        which = localRnd.nextInt(numPart);

        xTrial = correctPosition(Xs[which] + x);
        yTrial = correctPosition(Ys[which] + y);
        zTrial = correctPosition(Zs[which] + z);

        // Calculating the potential shift
        double oldE = 0, newE = 0;

        // firstly, old Energy
        for (int i = 0; i < numPart; i++) {
            oldE = oldE + getPotential(i, which);
        }

        // then, new Energy
        for (int i = 0; i < numPart; i++) {
            newE = newE + getPotential(i, which, xTrial, yTrial, zTrial);
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

    protected double getBoxSize() {
        return boxSize;
    }

    @Override
    /**
     * Used only by external watchers! gives back averaged potential (by last 5 configurations).
     */
    public double getEnergy() {
        return avgEnergy;
    }


    @Override
    public final boolean isFinished() {
        return finished;
    }


    @Override
    public int getNumPart() {
        return numPart;
    }

    @Override
    public int getT() {
        return T;
    }

    @Override
    public String getFolder() {
        return myFolder;
    }

    @Override
    public int getNumSteps() {
        return numSteps;
    }

    @Override
    public final void stop() {
        finished = true;
    }

    public String toString() {
        return "ensmbl:" + myFolder;
    }
}
