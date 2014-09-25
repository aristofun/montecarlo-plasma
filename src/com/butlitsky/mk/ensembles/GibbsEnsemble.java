package com.butlitsky.mk.ensembles;

import com.butlitsky.mk.options.CLOptions;
import com.butlitsky.mk.options.EOptions;
import org.apache.commons.math3.util.FastMath;

import java.io.IOException;
import java.util.ArrayDeque;
import java.util.Deque;

/**
 * Gibbs ensemble for e-i classical plasma.
 * <p/>
 * Straightforward implementation of ATHANASSIOS Z. PANAGIOTOPOULOS basic idea.
 * http://kea.princeton.edu/papers/varenna94/varenna.pdf
 * <p/>
 * Date: 23.09.14
 * Time: 18:20
 * <p/>
 * TODO: move all code duplication to separate NVTAssistant class to share with NVTEnsemble
 */
public abstract class GibbsEnsemble extends MetropolisEnsemble {

    // ----------–---------- File I/O tools ----------------------
    private GibbsConfigurationManager config;

    // ------------------- Common physics ------------------------
    private final int N; // total particles pairs number N = Nei*2
    private final int Nei; // total particles pairs number N = Nei*2
    private final double V; // total Volume (in Bohrs^3)
    protected final int T; // temperature (in K)
    private final double Gamma; // total initial point gamma

    // ---------–––––––– Common Monte-Carlo ----------------------

    // test particle movement coordinates
    private int trialType, trialIndex;
    private double xTrial, yTrial, zTrial;

    // current "accepted" state of the particles of both boxes
    private final double[][][] prtcls;

    // trial state for V1,V2 change step. To be switched with 'parts' on acceptance
    private final double[][][] tryPrtcls;


    // ------------ BOX mathematics ------------------------------

    // current box border index in the array of particles' positions
    // [0...boxborder) – BOX 1; [boxBorder...Nei) – BOX 2
    private int boxBorder;

    private double V1, V2; // in Bohrs^3
    private final double[] boxSize; // in Bohrs
    private final double[] halfBox; // in Bohrs
    private final int[] lengths; // [0] == boxBorder, [1] == Nei-boxBorder

    // must be recalclated on every volume chane
    private final double[] deltaX;
    private final double maxDXfactor;

    private float acceptance = 0; // acceptance rate
    private int acceptCnt = 0; // acceptance rate counter
    private int acceptIterations = 0;

    // --------------- Results accumulators -------------------------
    private final int avgPoints;
    private final double[] avgEnergies;

    private final double[] currentEnergies;

    // averaging energies values stack
    private final Deque<Double>[] energies;
    // density (1/Bohr^3) for past 'resolution' steps accumulator (each box)
    private final double[] densitiesSum;
    private final double[] densitiesAvg; // average densite for past 'resolution' steps

    private int densitiesIterations = 0;

    // ------------ Gibbs additional options ---------------------
    private final double maxDeltaV;

    // what move was made: 0. move random particle or 1. change V or 2. interchange particles
    private int lastStepType = 0;
    // what box was changed during last MC step: 0,1 or -1 (meaning both)
    private int lastBox = -1;

    protected GibbsEnsemble(EOptions options) {
        super(options,
              CLOptions.N_RESOLUTION_STEPS,
              CLOptions.N_RESOLUTION_STEPS * 3,
              CLOptions.N_RESOLUTION_STEPS * 5);

        N = opt.getNumParticles();
        Nei = opt.getNumParticles() / 2;
        T = opt.getT();
        V = N / (2 * opt.getDensity() * BOHR * BOHR * BOHR);

        Gamma = opt.getGamma();
        maxDXfactor = opt.getMaxDelta();
        maxDeltaV = CLOptions.MAX_DELTA_V;

        boxSize = new double[2];
        halfBox = new double[2];
        // we double 'cause total density is twice bigger
        V1 = V / 2;
        boxSize[0] = FastMath.cbrt(V1);

        V2 = V - V1;
        boxSize[1] = FastMath.cbrt(V2);

        halfBox[0] = boxSize[0] / 2.0;
        halfBox[1] = boxSize[1] / 2.0;

        // split particles initially by two parts
        boxBorder = Nei / 2;

        lengths = new int[2];
        lengths[0] = boxBorder;
        lengths[1] = Nei - boxBorder;

        // new in 8.0 – CLI params always multiplied by avgDistance
        deltaX = new double[2];
        deltaX[0] = (maxDXfactor == 0.0) ? boxSize[0] : maxDXfactor
                * FastMath.cbrt(1. / (2 * boxBorder / V1));
        deltaX[1] = (maxDXfactor == 0.0) ? boxSize[1] : maxDXfactor
                * FastMath.cbrt(1. / (2 * (Nei - boxBorder) / V2));

        System.out.println(
                myFolder + ": boxBorder=" + boxBorder + ", " +
                        "deltaV=" + SHORT_FORMAT.format(maxDeltaV));

        System.out.println(
                "gamma=" + SHORT_FORMAT.format(Gamma) + ", " +
                        "V(bohrs)=" + SHORT_FORMAT.format(V) + ", deltaX1=" + SHORT_FORMAT.format(
                        deltaX[0]));

        // initializing particles configuration: [e,i][X,Y,Z][0...Nei]
        prtcls = new double[2][3][Nei];
        tryPrtcls = new double[2][3][Nei];

        avgPoints = CLOptions.NUM_ENERGY_AVG_STEPS;
        avgEnergies = new double[2];

        densitiesSum = new double[2];
        densitiesAvg = new double[2];

        currentEnergies = new double[2];
        System.out.print(", AVG.=" + avgPoints);

        energies = new Deque[2];
        energies[0] = new ArrayDeque<>(avgPoints);
        energies[1] = new ArrayDeque<>(avgPoints);

        config = new GibbsConfigurationManager(this);
    }


    // --------------- State lifecycle management -----------------------------------

    /**
     * Must be called by children in the end of construction or by a factory.
     */
    public void loadState() {
        loadConfiguration();
        initEnergy();
        config.applyAdditionalStrategies();
    }

    private void loadConfiguration() {
        try {
            config.loadConfiguration();
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: can't find " + myFolder + " directory, worker is dead");
            stop();
            return;
        }
        System.out.println();
    }

    @Override
    protected void saveStateOnStop() {
        averageEnergies();

        config.saveConfiguration();
        config.saveLongTail();
        config.closeLongTail();
//        saveCorrelation();  todo: not fixed yet
    }

    private void averageEnergies() {
        for (int box = 0; box < 2; box++) {
            averageEnergy(box);
        }
    }

    // ------------------------------- Monte-Carlo game -------------------------------

    /**
     * performs one act of monte-karlo play randomly moves one particle and saves the state
     * with some probability
     */
    protected final boolean play(int step) {
        // 0. move random particle or 1. change V1 or 2. interchange particle
//        lastStepType = step % 3;
        lastStepType = 0;

        final double expoValue = testMove();

//        System.out.println("Expo val: " + expoValue);
        // transition probability checking
        // potential increased
        if (expoValue < 1) {        // potential increased, exp(-deltaE) < 1
            // compare the transition probability with random
            // All energies are in kT
            if (expoValue > myRandom(1.0)) {
                acceptTestMove();
                return true;
            }
        } else {     // potential decreased, accept configuration, exp(-deltaE) > 1
            acceptTestMove();
            return true;
        }

        return false;
    }

    /**
     * densities for past 'resolution' steps are accumulated to calculate averages
     */
    private void everyStepCalcs() {
        acceptIterations++;
        densitiesIterations++;

        densitiesSum[0] += boxBorder * 2 / V1;
        densitiesSum[1] += (Nei - boxBorder) * 2 / V2;
//        System.out.println("every step acceptIt/densIt: " + acceptIterations + " / " + densitiesIterations);
    }

    @Override
    protected void onTrialAccepted() {
        // what move was made: 0. move random particle or 1. change V or 2. interchange particles
        switch (lastStepType) {
            case 0:
                break;
//            case 0:
//                break;
//            case 0:
//                break;

            default:
                throw new IllegalStateException("Unknown MC test move type! Can't go on.");
        }

        if (lastBox > -1) {
            newEnergyStep(lastBox);
        } else {
            newEnergyStep(0);
            newEnergyStep(1);
        }

        acceptCnt++;
        everyStepCalcs();
    }

    @Override
    protected void onTrialRejected() {
        if (lastBox > -1) {
            oldEnergyStep(lastBox);
        } else {
            oldEnergyStep(0);
            oldEnergyStep(1);
        }

        acceptCnt--;
        everyStepCalcs();
    }

    @Override
    protected void doRareCalc() {
        System.out.println("i'm rare calc! step: " + getCurrStep());
        averageEnergies();
        config.saveConfiguration();
//        saveCorrelation(); todo: not fixed yet
    }

    @Override
    protected void doMidCalc() {
        System.out.println("i'm mid calc! step: " + getCurrStep());
//        config.calcCorrelation(); todo: not fixed yet
        config.saveLongTail();
    }

    @Override
    protected void doFrequentCalc() {
        // record acceptance rate & reset counter
        acceptance = ((1 + ((float) acceptCnt / acceptIterations)) / 2);
        acceptCnt = 0;
        System.out.println("frequent calc (acceptance/iterations): " + acceptance + " / " + acceptIterations);
        acceptIterations = 0;

        // record densities
        densitiesAvg[0] = densitiesSum[0] / densitiesIterations;
        densitiesAvg[1] = densitiesSum[1] / densitiesIterations;
        densitiesIterations = 0;
        densitiesSum[0] = 0.0;
        densitiesSum[1] = 0.0;

//        writeAdditionalStateFile();
    }

    private final void acceptTestMove() {
        // what move was made: 0. move random particle or 1. change V or 2. interchange particles
        switch (lastStepType) {
            case 0:
                acceptParticleMove();
                break;
//            case 1:
            default:

        }
    }

    /**
     * Actually updates all 'trial' particle coordinates in the given box
     */
    private final void acceptParticleMove() {
        prtcls[trialType][0][trialIndex] = xTrial;
        prtcls[trialType][1][trialIndex] = yTrial;
        prtcls[trialType][2][trialIndex] = zTrial;
    }


    // ------------- Averages calculation (energies, densities etc.) --------------------

    /**
     * sets initial average energy for current configuration
     */
    private void initEnergy() {
        for (int i = 0; i < 2; i++) {
            if (avgEnergies[i] == 0) {
                newEnergyStep(i);
                averageEnergy(i);
            } else {
                currentEnergies[i] = avgEnergies[i];
            }
        }

    }

    /**
     * Record current energy value to averager array
     */

    private final void newEnergyStep(int box) {
        currentEnergies[box] = getCurrentEnergy(box);
        oldEnergyStep(box);
    }

    private final void oldEnergyStep(int box) {
        if (energies[box].size() > avgPoints - 1) energies[box].pollLast();
        energies[box].addFirst(currentEnergies[box]);
    }

    private final void averageEnergy(int box) {
        double enrg = 0.0;
        for (Double val : energies[box]) {
            enrg = enrg + val;
        }
        avgEnergies[box] = enrg / energies[box].size();
    }

    /**
     * Calculates and returns full system energy from scratch, based on current particles
     * configuration. No changes made to configuration.
     *
     * @param whichBox (0 for first box, 1 - for second)
     */
    protected final double getCurrentEnergy(final int whichBox) {
        double newEn = 0;
        final int length = lengths[whichBox];
        final int offset = (whichBox == 0) ? 0 : boxBorder;

        for (int i = 0; i < length * 2; i++) {
            for (int j = i + 1; j < length * 2; j++) {
                newEn = newEn +
                        getEnergy(
                                i / length, (i % length) + offset,
                                j / length, (j % length) + offset,
                                whichBox);
            }
        }
        return newEn;
    }


    /**
     * Make one of three types of trial moves.
     *
     * @return expo value for MC probability weighting
     */
    private final double testMove() {
        // what move was made: 0. move random particle or 1. change V or 2. interchange particles
        switch (lastStepType) {
            case 0:
                //                System.out.println("delta: " + delta);
                return FastMath.exp(-testMoveParticle());
            case 1:
//                return testChangeV();
            case 2:
//                return testSwitchParticles();
            default:
                throw new IllegalStateException("Unknown MC test move type! Can't go on.");
        }
    }
    // --------------- PARTICELS movement routines, boundary conditions -------------------

    /**
     * returns trial potential shift for moved particle
     */
    private final double testMoveParticle() {
        // general parameters first set
        lastBox = nextInt(2);
        final int length = lengths[lastBox];
        final int which = nextInt(length * 2);
        final int offset = (lastBox == 0) ? 0 : boxBorder;

        // setting new trial index and coordinates
        trialType = which / length;
        trialIndex = (which % length) + offset;

        final double x = myRandom() * deltaX[lastBox];
        final double y = myRandom() * deltaX[lastBox];
        final double z = myRandom() * deltaX[lastBox];

        xTrial = correctPosition(prtcls[trialType][0][trialIndex] + x, boxSize[lastBox]);
        yTrial = correctPosition(prtcls[trialType][1][trialIndex] + y, boxSize[lastBox]);
        zTrial = correctPosition(prtcls[trialType][2][trialIndex] + z, boxSize[lastBox]);

        // Calculating the potential shift
        double oldE = 0, newE = 0;

        // firstly, old Energy
        oldE = sumPotential(trialType, trialIndex,
                            prtcls[trialType][0][trialIndex],
                            prtcls[trialType][1][trialIndex],
                            prtcls[trialType][2][trialIndex],
                            prtcls,
                            lastBox,
                            length, offset
        );

        // then, new Energy
        newE = sumPotential(trialType, trialIndex, xTrial, yTrial, zTrial, prtcls, lastBox,
                            length, offset);


        return newE - oldE;
    }

    /**
     * Sum of the potentials between given particle (with given coords) and all other particles in
     * the box
     *
     * @param length
     * @param offset
     * @return
     */
    private double sumPotential(
            final int particleType, final int particleIndex,
            final double x, final double y, final double z,
            final double[][][] particls,
            final int whichBox,
            int length, int offset
    ) {
        double result = 0.0;

        for (int type = 0; type < 2; type++) {
            for (int i = 0; i < length; i++) {
                result = result + getPotential(type, i + offset, particleType, particleIndex,
                                               x, y, z, particls, whichBox);
            }
        }
        return result;
    }


    private final double dSquared(double dx, double dy, double dz, final int whichBox) {
        final double halfbox = halfBox[whichBox];
        dx = fit2box(dx, halfbox);
        dy = fit2box(dy, halfbox);
        dz = fit2box(dz, halfbox);
        return ((dx * dx) + (dy * dy) + (dz * dz));
    }

    private final double fit2box(double dx, final double halfbox) {
//        dx = abs(dx);
//        return (dx > halfBox) ? (halfBox - dx % halfBox) : dx;
        if (dx > halfbox) {
            dx -= (halfbox + halfbox);
        } else if (dx < -halfbox) {
            dx += (halfbox + halfbox);
        }
        return dx;
    }

    /**
     * set right position inside the Box
     */
    private final double correctPosition(final double coord, final double box_size) {
        return ((box_size + coord % box_size) % box_size);
    }


    // ----------------- Potential & Energies ----------------------------------

    /**
     * Two particle potential value between given particles with given coordinates.
     *
     * @param type1,2  – first particle type in prtcls[type][][] array
     * @param num1     – first particle number in prtcls[][][num] array
     * @param particls - which particles array to use (either main or test)
     * @param whichBox - which box parameters to use (boxSize, halfSize, etc.)
     * @return energy value in kT units
     */
    private final double getPotential(
            final int type1, final int num1,
            final int type2, final int num2,
            final double trialX, final double trialY, final double trialZ,
            final double[][][] particls,
            final int whichBox
    ) {

        if (type1 == type2 && num1 == num2) {
            return 0;
        }

        final double r = Math.sqrt(
                dSquared(
                        particls[type1][0][num1] - trialX,
                        particls[type1][1][num1] - trialY,
                        particls[type1][2][num1] - trialZ,
                        whichBox)
        );

        if (type1 == 0)   // Electrons
        {
            if (type2 == 0) // E-E
                return getPotentialAsym(r, true, false);
            else              // E-I
                return getPotentialAsym(r, false, false);
        } else                   // Ions
        {
            if (type2 == 0) // I-E
                return getPotentialAsym(r, false, false);
            else              // I-I
                return getPotentialAsym(r, false, true);
        }
    }

    protected abstract double getPotentialAsym(double r, boolean ee, boolean ii);

    /**
     * Two particle energy value between particles with given coordinates.
     *
     * @param type1,2  – first particle type in prtcls[type][][] array
     * @param num1,2   – first particle number in prtcls[][][num] array
     * @param whichBox - which box parameters to use (boxSize, halfSize, etc.)
     * @return energy value in kT units
     */
    private final double getEnergy(final int type1, final int num1,
                                   final int type2, final int num2,
                                   final int whichBox) {

        final double r = Math.sqrt(
                dSquared(
                        prtcls[type1][0][num1] - prtcls[type2][0][num2],
                        prtcls[type1][1][num1] - prtcls[type2][1][num2],
                        prtcls[type1][2][num1] - prtcls[type2][2][num2],
                        whichBox)
        );

        if (type1 == type2 && num1 == num2) {
            return 0;
        }

        if (type1 == 0)   // Electrons
        {
            if (type2 == 0) // E-E
                return getEnergyAsym(r, true, false);
            else              // E-I
                return getEnergyAsym(r, false, false);
        } else                   // Ions
        {
            if (type2 == 0) // I-E
                return getEnergyAsym(r, false, false);
            else              // I-I
                return getEnergyAsym(r, false, true);
        }
    }

    protected abstract double getEnergyAsym(double r, boolean ee, boolean ii);


    // -------------- Public contracts ------------------------------------------------

    /**
     * @return gamma1 E1 gamma2 E2 values for plotting two separate boxes results
     */
    //        gamma = e * e * FastMath.cbrt(2 * density) / (k * T);
    @Override
    public double[] getCurrentResult() {
        // get 'n' in 1/cm^3
        final double n1 = densitiesAvg[0] / (BOHR * BOHR * BOHR);
        final double n2 = densitiesAvg[1] / (BOHR * BOHR * BOHR);
        // get gammas
        final double gamma1 = e * e * FastMath.cbrt(n1) / (k * T);
        final double gamma2 = e * e * FastMath.cbrt(n2) / (k * T);

        return new double[]{
                gamma1,
                avgEnergies[0] / N,
                gamma2,
                avgEnergies[1] / N,
        };
    }

    @Override
    public int getT() { return T; }

    // --------- SETTERS for external configuration loading & writing ----------------
    void setBoxBorder(int border) {
        boxBorder = border;
    }

    int getBoxBorder() {
        return boxBorder;
    }

    double[] getBoxSizes() {return boxSize; }

    void setAvgEnergies(double en1, double en2) {
        avgEnergies[0] = en1;
        avgEnergies[1] = en2;
    }

    double getAvgEnergy1() {return avgEnergies[0];}

    double getAvgEnergy2() {return avgEnergies[1];}

    double[][][] getParticles() {return prtcls;}

}
