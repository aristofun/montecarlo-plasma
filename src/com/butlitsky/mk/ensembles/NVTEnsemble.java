package com.butlitsky.mk.ensembles;

import com.butlitsky.mk.options.CLOptions;
import com.butlitsky.mk.options.EOptions;
import org.apache.commons.math3.util.FastMath;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Deque;
import java.util.List;

import static org.apache.commons.math3.util.FastMath.ceil;
import static org.apache.commons.math3.util.FastMath.exp;

/**
 * Basic class for all the metropolis abstractions. Assuming NVT single box ensemble.
 * <p/>
 * Date: 03.03.13
 * Time: 13:04
 */
public abstract class NVTEnsemble extends MetropolisEnsemble {

    private final int avgPoints;

    private Path myConfigPath;
    private Path myCorrPath;

    private final boolean saveLongTail;
    private BufferedWriter longTailWriter;


    private final int numPart;
    private int which;
    private double xTrial;
    private double yTrial;
    private double zTrial;

    // ------------ mathematics --------------
    protected final int T;
    private final double boxSize;
    private final double halfBox;
    private final double maxDelta;

    // -–––– correlation stuff ---------------
    private final double corrNormirovka;
    private final double corrDr;
    private final double[][] corrArray = new double[3][];
    private int corrAverager = 0;

    protected final double[] Xs;
    protected final double[] Ys;
    protected final double[] Zs;

    // averaging energies values stack
    private final Deque<Double> energies;
    private double avgEnergy = 0;
    private double currentEnergy = 0;


    protected NVTEnsemble(EOptions options) {
        super(options,
              options.getNumParticles() + 1,
              options.getNumParticles() * 3 + 7,
              options.getNumParticles() * 7 + 11);

        numPart = opt.getNumParticles();
        T = opt.getT();

        boxSize = FastMath.cbrt(numPart / (2 * opt.getDensity())) / BOHR; //
        // parameter is Ne(Ni),
        // we double 'cause total density is twice bigger

        halfBox = boxSize / 2.0;
        // average e-i distance => x 2
        double avgDistance = FastMath.cbrt(1. / (2. * opt.getDensity())) / BOHR;
        double factor = opt.getMaxDelta();

        // new in 8.0 – CLI params always multiplied by avgDistance
        maxDelta = (factor == 0.0) ? boxSize : factor * avgDistance;

        System.out.print(myFolder + ": avg.dist=" + SHORT_FORMAT.format(avgDistance) + ", " +
                                 "delta=" + SHORT_FORMAT.format(maxDelta) + ", boxSize="
                                 + SHORT_FORMAT.format(boxSize));

        corrNormirovka = 4. * (boxSize * boxSize * boxSize) / (numPart * numPart);
        corrDr = Math.sqrt(3.0 * boxSize * boxSize) / 1.99999999999 / (CORR_LENGTH);

        corrArray[0] = new double[CORR_LENGTH];
        corrArray[1] = new double[CORR_LENGTH];
        corrArray[2] = new double[CORR_LENGTH];

        // initiating particles configuration
        Xs = new double[numPart];
        Ys = new double[numPart];
        Zs = new double[numPart];

        // OPTIONS first bit == save longtail
        saveLongTail = ((options.getStrategy() & 1) == 1);
        avgPoints = (CLOptions.NUM_ENERGY_AVG_STEPS < 0) ?
                getNumSteps() - CLOptions.INITIAL_STEPS : CLOptions.NUM_ENERGY_AVG_STEPS;
        System.out.print(", AVG.=" + avgPoints);

        energies = new ArrayDeque<>(avgPoints);
    }

    /**
     * Must be called by children in the end of construction. So Ensemble is ready to run right
     * after it's successfully created.
     */
    public void loadState() {
        loadFromStateFile();
        initEnergy();
        applyAdditionalStrategies();

        if (!opt.isOld()) {      // save initially creted config for the first time
            saveConfiguration();
        }
    }

    /**
     * sets average potential from file or resets it to zero
     */
    private void initEnergy() {

        if (avgEnergy == 0) {
            newEnergyStep();
            averageEnergy();
        } else {
            currentEnergy = avgEnergy;
        }
    }


    private void applyAdditionalStrategies() {
        if (saveLongTail) {
            try {
                longTailWriter = Files.newBufferedWriter(GibbsConfigurationManager.getPath(myFolder + "/" + LONGTAIL_FILE),
                                                         Charset.forName("UTF-8"),
                                                         StandardOpenOption.CREATE,
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
        myConfigPath = GibbsConfigurationManager.getPath(myFolder + "/" + STATE_FILE);
        myCorrPath = GibbsConfigurationManager.getPath(myFolder + "/" + CORR_FILE);

        try {
            Files.createDirectories(GibbsConfigurationManager.getPath(myFolder));
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: can't find " + myFolder + " directory, worker is dead");
            stop();
            return;
        }

        // reading old state
        if (opt.isOld()) {
            if (Files.exists(myConfigPath)) {
                try {
                    loadArrays(Files.readAllLines(myConfigPath, Charset.forName("UTF-8")));
                } catch (Exception e) {
                    System.out.println("WARNING: failed to read config for " + myFolder);
                    System.out.println(e.getLocalizedMessage());
                    opt.setOld(false);
                }
            } else {
                opt.setOld(false);
            }
        }

        if (!opt.isOld()) {
            System.out.print("! from scratch");
            initParticlesPosition(); // initial distribution, reset counters
        }

        System.out.println();
    }

    /**
     * fill random arrays of particles, resetting all counters
     */
    private void initParticlesPosition() {
        // uniformly distribute particles in the box NaCL structure (fcc)
        if (CLOptions.START_FROM_FCC) {
            final int Nx = (int) ceil(Math.cbrt(numPart));
            final double Lx = boxSize / Nx;

            for (int j = 0; j < numPart; j++) {
                int xLayer = j / (Nx * Nx);
                int yLayer = (j % (Nx * Nx)) / Nx;
                int zLayer = (j % (Nx * Nx)) % Nx;

                if ((j % 2) == 0) { // ion, electron, ion, electron...
                    Xs[j / 2] = xLayer * Lx;
                    Ys[j / 2] = yLayer * Lx;
                    Zs[j / 2] = zLayer * Lx;
                } else {
                    Xs[numPart / 2 + j / 2] = xLayer * Lx;
                    Ys[numPart / 2 + j / 2] = yLayer * Lx;
                    Zs[numPart / 2 + j / 2] = zLayer * Lx;
                }

            }
        } else {
            // random configuration.
            for (int j = 0; j < numPart; j++) {
            /* spreading the particles */
                Xs[j] = myRandom(boxSize);
                Ys[j] = myRandom(boxSize);
                Zs[j] = myRandom(boxSize);
            }
        }
    }

    /**
     * Record current energy value to averager array
     */
    private final void newEnergyStep() {
        currentEnergy = getCurrentEnergy();
        oldEnergyStep();
    }

    private final void oldEnergyStep() {
        if (energies.size() > avgPoints - 1) energies.pollLast();
        energies.addFirst(currentEnergy);
    }

    private final void averageEnergy() {
        double enrg = 0.0;

        for (Double val : energies) {
            enrg = enrg + val;
        }

        avgEnergy = enrg / energies.size();
    }

    /**
     * calculates full system potential from scratch, based on current particles configuration
     */
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
        final double r = Math.sqrt(dSquared(trialX - Xs[i], trialY - Ys[i], trialZ - Zs[i]));

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

    protected abstract double getPotentialAsym(double r, boolean ee, boolean ii);

    private final double getEnergy(int i, int j) {
        final double r = Math.sqrt(dSquared(Xs[j] - Xs[i], Ys[j] - Ys[i], Zs[j] - Zs[i]));

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

    protected abstract double getEnergyAsym(double r, boolean ee, boolean ii);

    private final double dSquared(double dx, double dy, double dz) {
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
     * throws IOException if no valid config could be read
     */
    private void loadArrays(List<String> strings) throws Exception {
        // first line: curr_step potential
        String step_en = strings.remove(0);
        // #step, average potential, ...

        setCurrStep(Integer.parseInt(step_en.split("\\s+")[0]));
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

    @Override
    protected void saveStateOnStop() {
        saveConfiguration();
        saveCorrelation();
        saveLongTail();
        closeLongTail();
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
            Files.write(myCorrPath, strings, Charset.forName("UTF-8"));
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: failed to save correlation for " + myFolder);
        }
    }

    private void saveConfiguration() {
        averageEnergy();

        // create strings list from arrays
        try {
            BufferedWriter writer = Files.newBufferedWriter(myConfigPath, Charset.forName("UTF-8"));

            writer.write("" + getCurrStep() + "\t"
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


    private void closeLongTail() {
        if (saveLongTail) try {
            longTailWriter.flush();
            longTailWriter.close();
        } catch (IOException e1) {
            System.out.println("ERROR: can't close long tail writer for " + myFolder);
        }
    }

    private void saveLongTail() {
        if (saveLongTail) try {
            writeStateTo(longTailWriter);
            longTailWriter.flush();
        } catch (IOException e1) {
            System.out.println("ERROR: can't write " + myFolder + " long tail to writer!");
        }

    }

    private final void fillCorrAray(int i, int j, int corrIndex) {
        if (i != j) {
            // get the index of a radius in array
            final int idx =
                    (int) (Math.sqrt(dSquared(Xs[i] - Xs[j], Ys[i] - Ys[j],
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
     * @param step
     */
    protected final boolean play(int step) {
        // move random particle
        final double deltaE = moveParticle();
        // transition probability checking
//        System.out.println("deltaE/exp(): " + deltaE + " / " + exp(-deltaE));
        // potential increased
        if (deltaE > 0) {
            // compare the transition probability with random
            // All energies are in kT
            if (exp(-deltaE) >= myRandom(1.00000001)) {
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

    @Override
    protected void onTrialRejected() {
        oldEnergyStep();
    }

    @Override
    protected void onTrialAccepted() {
        newEnergyStep();
    }

    @Override
    protected void doRareCalc() {
        averageEnergy();
        saveConfiguration();
        saveCorrelation();
    }

    @Override
    protected void doMidCalc() {
        calcCorrelation();
        saveLongTail();
    }

    @Override
    protected void doFrequentCalc(int curr_step) {
        // do nothing
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
        which = nextInt(numPart);

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
    }

    protected double getBoxSize() {
        return boxSize;
    }

    @Override
    public double[] getCurrentResult() {
        return new double[]{avgEnergy / numPart};
    }

    public int getNumPart() {
        return numPart;
    }

    @Override
    public int getT() {
        return T;
    }
}
