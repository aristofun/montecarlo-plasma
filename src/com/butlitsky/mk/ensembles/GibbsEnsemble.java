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
    public static final String GIBBS_STATE_FILE = "gibbs2box.dat";


    // ----------–---------- File I/O tools ----------------------
    private GibbsConfigurationManager config;

    // ------------------- Common physics ------------------------
    private final int N; // total particles number N = Nei*2
    private final int Nei; // total particles pairs number N = Nei*2
    private final double Volume; // total Volume (in cm^-3)
    protected final int T; // temperature (in K)
    private final double Gamma; // total initial point gamma

    // ---------–––––––– Common Monte-Carlo ----------------------

    // test particle movement coordinates
    private int trialType, trialIndex;
    // trial coordinates for pair of particles (in case of particles switch trial movement
    private final double[] xTrial = new double[2], yTrial = new double[2], zTrial = new double[2];
    private final double[] trialBoxSize = new double[2];
    // trial Volume 1 new value (in cm-3)
    private double trialV0;

    // current "accepted" state of the particles of both boxes
    private double[][][] prtcls;

    // trial state for V1,V2 change step. To be switched with 'parts' on acceptance
    private double[][][] trialPrtcls;

    private double[][][] tmp_pointer;


    // ------------ BOX mathematics ------------------------------

    // current box border index in the array of particles' positions
    // [0...boxborder) – BOX 1; [boxBorder...Nei) – BOX 2
    private int boxBorder;

    private final double[] V = new double[2]; // V1, V2 in cm^-3
    private final double[] boxSize = new double[2]; // in Bohrs
    private final double[] halfBox = new double[2]; // in Bohrs
    private final int[] lengths = new int[2]; // [0] == boxBorder, [1] == Nei-boxBorder

    // must be recalclated on every volume chane
    private final double[] deltaX = new double[2];
    private final double maxDXfactor;

    private float acceptance = 0; // acceptance rate
    private int acceptCnt = 0; // acceptance rate counter
    private int acceptIterations = 0;

    // --------------- Results accumulators -------------------------
    private final int avgPoints;
    private final double[] reducedEnrgyAvg = new double[2];
    private final Deque<Double>[] reducedEnergies = new Deque[2];

    private final double[] currentEnergy = new double[2];

    // density (1/cm^-3) for past 'resolution' steps accumulator (each box)
    private final double[] densitiesSum = new double[2];
    private final double[] densitiesAvg = new double[2]; // average densite for past 'resolution' steps
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
        Volume = N / (2 * opt.getDensity());       // V is in CM^3 !!!

        Gamma = opt.getGamma();
        maxDXfactor = opt.getMaxDelta();
        maxDeltaV = CLOptions.MAX_DELTA_V;

        // we double 'cause total density is twice bigger
        V[0] = Volume / 2.;
        V[1] = Volume - V[0];

        updateBoxSizes(false);

        // split particles initially by two parts
        boxBorder = Nei / 2;

        updateLengths();
        updateDeltaX();


        System.out.println(myFolder + ": boxBorder=" + boxBorder + ", " +
                                   "deltaV=" + SHORT_FORMAT.format(maxDeltaV));

        System.out.println("gamma=" + SHORT_FORMAT.format(Gamma) + ", " + "V(cm-3)=" +
                                   SHORT_FORMAT.format(Volume) + ", deltaX1="
                                   + SHORT_FORMAT.format(deltaX[0]));

        // initializing particles configuration: [e,i][X,Y,Z][0...Nei]
        prtcls = new double[2][3][Nei];
        trialPrtcls = new double[2][3][Nei];

        avgPoints = CLOptions.NUM_ENERGY_AVG_STEPS;
        System.out.print(", AVG.=" + avgPoints);

        reducedEnergies[0] = new ArrayDeque<>(avgPoints);
        reducedEnergies[1] = new ArrayDeque<>(avgPoints);

        config = new GibbsConfigurationManager(this);
    }

    /**
     * Call only after V[] changed!
     *
     * @param afterTestV if true – trialBoxSize[] is used as a source,
     *                   otherwise - calculated based on current V[] data
     */
    private final void updateBoxSizes(boolean afterTestV) {
        if (afterTestV) {
            boxSize[0] = trialBoxSize[0];
            boxSize[1] = trialBoxSize[1];
        } else {
            boxSize[0] = FastMath.cbrt(V[0]) / BOHR;
            boxSize[1] = FastMath.cbrt(V[1]) / BOHR;
        }
        halfBox[0] = boxSize[0] / 2.0;
        halfBox[1] = boxSize[1] / 2.0;
    }

    /**
     * Call only after {@link #updateBoxSizes(boolean)} and {@link #updateLengths()}!
     */
    private final void updateDeltaX() {// new in 8.0 – CLI params always multiplied by avgDistance
        deltaX[0] = (maxDXfactor == 0.0) ? boxSize[0] : maxDXfactor
                * FastMath.cbrt(V[0] / (2 * lengths[0])) / BOHR;

        deltaX[1] = (maxDXfactor == 0.0) ? boxSize[1] : maxDXfactor
                * FastMath.cbrt(V[1] / (2 * lengths[1])) / BOHR;
    }

    /**
     * Call only after boxBorder updated!
     */
    private final void updateLengths() {
        lengths[0] = boxBorder;
        lengths[1] = Nei - boxBorder;
    }


    // --------------- State lifecycle management -----------------------------------

    /**
     * Must be called by children in the end of construction or by a factory.
     */
    public void loadState() {
        loadConfiguration();

        updateLengths();
        updateDeltaX();

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
        config.workOnMidCalc();
        config.closeLongTails();
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
        lastStepType = step % 3;

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
        //  density is in cm^-3 !!!
        densitiesSum[0] += lengths[0] * 2 / V[0];
        densitiesSum[1] += lengths[1] * 2 / V[1];
    }

    @Override
    protected void onTrialAccepted() {
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
        config.saveConfiguration();
//        saveCorrelation(); todo: not fixed yet
    }

    @Override
    protected void doMidCalc() {
        averageEnergies();
        System.out.println("i'm mid calc! step: " + getCurrStep());
//        config.calcCorrelation(); todo: not fixed yet
        config.workOnMidCalc();
    }

    @Override
    protected void doFrequentCalc(int curr_step) {
        // record acceptance rate & reset counter
        acceptance = ((1 + ((float) acceptCnt / acceptIterations)) / 2);
        acceptCnt = 0;
        System.out.println("frequent calc (acceptance/iterations): " + acceptance + " / " + acceptIterations);
        acceptIterations = 0;

        // record densities
        densitiesAvg[0] = densitiesSum[0] / densitiesIterations;
        densitiesAvg[1] = densitiesSum[1] / densitiesIterations;
        densitiesIterations = 1;
        densitiesSum[0] = densitiesAvg[0];
        densitiesSum[1] = densitiesAvg[1];

        config.workOnFrequentCalc(curr_step);
    }

    private final void acceptTestMove() {
        // what move was made: 0. move random particle or 1. change V or 2. interchange particles
        switch (lastStepType) {
            case 0:
                acceptParticleMove();
                break;
            case 1:
                acceptVolumeChange();
                break;
            case 2:
                acceptParticleSwitch(); // must update length[], deltaX[] and boxBorder and
                // actually switch particles between boxes
                break;
            default:
                throw new IllegalStateException("trial move type not supported");
        }
    }

    /**
     * deltaX[], volumes, boxSizes and densities updated,
     * as well as prtcls array switched to new trialPrtcls.
     */
    private void acceptVolumeChange() {
        tmp_pointer = prtcls;
        prtcls = trialPrtcls;
        trialPrtcls = tmp_pointer;
        tmp_pointer = null;

        V[0] = trialV0;
        V[1] = Volume - trialV0;

        updateBoxSizes(true);
        updateDeltaX();
    }

    /**
     * Updates length[], deltaX[], boxBorder, densities and modifies prtcls array.
     */
    private void acceptParticleSwitch() {
//        trialIndex; // which pair
//        xTrial, yTrial, zTrial; // their new coords

        // [0... trialIndex ... ][boxBorder ... Nei]
        // =>
        // [0...  ...] [boxBorder-1 ... trialpartlc ... Nei]

        final int src_pos, dest_pos, length, deltaBB;

        if (trialIndex < boxBorder) {
            src_pos = trialIndex + 1;
            dest_pos = trialIndex;
            length = boxBorder - trialIndex - 1;
            deltaBB = -1;
        } else {
            src_pos = boxBorder;
            dest_pos = boxBorder + 1;
            length = trialIndex - boxBorder;
            deltaBB = 0;
        }

//        System.out.println("indx: " + trialIndex);
//        System.out.println("before: " + Arrays.toString(prtcls[0][0]));
        System.arraycopy(prtcls[0][0], src_pos, prtcls[0][0], dest_pos, length);
        System.arraycopy(prtcls[0][1], src_pos, prtcls[0][1], dest_pos, length);
        System.arraycopy(prtcls[0][2], src_pos, prtcls[0][2], dest_pos, length);
        System.arraycopy(prtcls[1][0], src_pos, prtcls[1][0], dest_pos, length);
        System.arraycopy(prtcls[1][1], src_pos, prtcls[1][1], dest_pos, length);
        System.arraycopy(prtcls[1][2], src_pos, prtcls[1][2], dest_pos, length);

        prtcls[0][0][boxBorder + deltaBB] = xTrial[0];
        prtcls[0][1][boxBorder + deltaBB] = yTrial[0];
        prtcls[0][2][boxBorder + deltaBB] = zTrial[0];
        prtcls[1][0][boxBorder + deltaBB] = xTrial[1];
        prtcls[1][1][boxBorder + deltaBB] = yTrial[1];
        prtcls[1][2][boxBorder + deltaBB] = zTrial[1];

        boxBorder = (deltaBB == 0) ? boxBorder + 1 : boxBorder - 1;
        updateLengths();
        updateDeltaX();
    }

    /**
     * Actually updates all 'trial' particle coordinates in the given box
     */
    private final void acceptParticleMove() {
        prtcls[trialType][0][trialIndex] = xTrial[0];
        prtcls[trialType][1][trialIndex] = yTrial[0];
        prtcls[trialType][2][trialIndex] = zTrial[0];
    }


    // ------------- Averages calculation (energies, densities etc.) --------------------

    /**
     * sets initial average energy for current configuration
     */
    private void initEnergy() {
        for (int i = 0; i < 2; i++) {
            if (currentEnergy[i] == 0) {
                newEnergyStep(i);
            } else {
                oldEnergyStep(i);
            }
            averageEnergy(i);
        }
    }

    /**
     * Record current energy value to averager array
     */
    private final void newEnergyStep(int box) {
        currentEnergy[box] = getCurrentEnergy(box) / (2. * lengths[box]);
        oldEnergyStep(box);
    }

    private final void oldEnergyStep(int box) {
        if (reducedEnergies[box].size() > avgPoints - 1) reducedEnergies[box].pollLast();
        reducedEnergies[box].addFirst(currentEnergy[box]);
    }


    /**
     * Total potential energy of a box inside given particls array.
     * Global boxBorder and lengths[] are used (assuming prtcls and trialPrtcls are different in
     * volumes only).
     *
     * @param whichBox
     * @return potential energy for current coords
     */
    private final double getCurrentPotential(final double[][][] particls, final int whichBox,
                                             final double half_box) {
        double newPot = 0.0;
        final int length = lengths[whichBox];
        final int offset = whichBox * boxBorder;


        for (int i = 0; i < length * 2; i++) {
            for (int j = i + 1; j < length * 2; j++) {
                final int x1 = i / length;
                final int x2 = j / length;
                final int y1 = (i % length) + offset;
                final int y2 = (j % length) + offset;

                newPot = newPot +
                        getPotential(x1,
                                     particls[x1][0][y1],
                                     particls[x1][1][y1],
                                     particls[x1][2][y1],
                                     x2,
                                     particls[x2][2][y2],
                                     particls[x2][2][y2],
                                     particls[x2][2][y2],
                                     half_box
                        );
            }
        }
        return newPot;
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


    private final void averageEnergy(int box) {
        double enrg = 0.0;
        for (Double val : reducedEnergies[box]) {
            enrg = enrg + val;
        }
        reducedEnrgyAvg[box] = enrg / reducedEnergies[box].size();
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
                return FastMath.exp(-testMoveParticle());
            case 1:
                return testChangeV();
            case 2:
                return testSwitchParticles();
            default:
                throw new IllegalStateException("Unknown MC test move type! Can't go on.");
        }
    }

    // --------------- PARTICELS movement routines, boundary conditions -------------------


    /**
     * deltaX[], volumes, boxSizes and densities must be updated on accepted move,
     * as well as prtcls array switched to new trialPrtcls.
     *
     * @return probability coefficient for trial volume change of the Volume 1
     */
    private final double testChangeV() {
        lastBox = -1;
        final double deltaV = myRandom() * maxDeltaV * FastMath.min(V[0], V[1]);
//        final double deltaV = 0.2 * maxDeltaV * FastMath.min(V[0], V[1]); xxx debug hack
        trialV0 = V[0] + deltaV;
        trialBoxSize[0] = FastMath.cbrt(trialV0) / BOHR;
        trialBoxSize[1] = FastMath.cbrt(Volume - trialV0) / BOHR;

        // 0. fill trialPrtcls array
        scaleTrialPrtcls();

        // 1. find prtcls and trialPrtcls array potential energies
        final double prtclsE = getCurrentPotential(prtcls, 0, halfBox[0])
                + getCurrentPotential(prtcls, 1, halfBox[1]);

        final double trialE = getCurrentPotential(trialPrtcls, 0, trialBoxSize[0] / 2.0)
                + getCurrentPotential(trialPrtcls, 1, trialBoxSize[1] / 2.0);

        final double NVcoeff =
                Math.log(trialV0 / V[0]) * lengths[0] * 2
                        + Math.log((V[1] - deltaV) / V[1]) * lengths[1] * 2;

        double res = FastMath.exp(NVcoeff + prtclsE - trialE);

        // XXX todo: remove after debug
//        System.out.println("deltaE:  " + (trialE - prtclsE) + ", deltaV: " + deltaV);
//        System.out.println("volumes:  " + V[0]/V[1]);
//        System.out.println("NVcoeff: " + NVcoeff + ", deltaV: "
//                                   + (int) (100 * deltaV / FastMath.min(V[0], V[1]))
//        + ", [ " + lengths[0] + " / " + lengths[1] + " ]");
//        System.out.println("Exp:     " + res);

//        try {
//            Thread.sleep(100);
//        } catch (InterruptedException e1) {
//            e1.printStackTrace();
//        }
        return res;
    }

    /**
     * Copies all prtcs coords to trialPrtcls with scaled coords, according to trialBoxSize[] data
     */
    private void scaleTrialPrtcls() {
        final double[] scalingCoeff = new double[]{
                trialBoxSize[0] / boxSize[0],
                trialBoxSize[1] / boxSize[1]};

        // scale all the particles coords (copy scaled coordinates to trialPrtcls array)
        for (int type = 0; type < 2; type++) {
            for (int i = 0; i < Nei; i++) {
                copyCoords(scalingCoeff[(i < boxBorder) ? 0 : 1], type, i);
            }
        }
    }

    /**
     * Copy given XYZ coords from prtcls -> trialPrtcls with given scaling factor
     *
     * @param scalingCoeff
     * @param type
     * @param i
     */
    private final void copyCoords(final double scalingCoeff, final int type, final int i) {
        trialPrtcls[type][0][i] = prtcls[type][0][i] * scalingCoeff;
        trialPrtcls[type][1][i] = prtcls[type][1][i] * scalingCoeff;
        trialPrtcls[type][2][i] = prtcls[type][2][i] * scalingCoeff;
    }

    /**
     * trialIndex is used as a random particle pair pointer.
     * lengths[], deltaX[], boxBorder and densities must be updated on accepted move.
     *
     * @return probability coefficient for single particles pair switch between box 1 and box 2
     */
    private final double testSwitchParticles() {
        lastBox = -1;   // both boxes are touched
        trialIndex = nextInt(Nei); // which particles pair should be moved
        final int fromBox = (trialIndex < boxBorder) ? 0 : 1;
        final int toBox = 1 - fromBox;

        // (N2 - 2)*V1/(N1 + 2)*V2 – transfer from 2 -> 1
        final double expoCoefficient =
                (lengths[fromBox] - 1) * V[toBox] / (V[fromBox] * (lengths[toBox] + 1));

        // XXX todo: remove after debug
//        long start = System.currentTimeMillis();
//        System.out.println("trial/boxBorder: " + trialIndex + "/" + boxBorder
//                                   + ", expoCoeff: " + expoCoefficient);

        // new random test coordinates of particle pair in new box
        xTrial[0] = myRandom(boxSize[toBox]);
        yTrial[0] = myRandom(boxSize[toBox]);
        zTrial[0] = myRandom(boxSize[toBox]);
        xTrial[1] = myRandom(boxSize[toBox]);
        yTrial[1] = myRandom(boxSize[toBox]);
        zTrial[1] = myRandom(boxSize[toBox]);


        // delta energy of particle pair transfer - !!! with inverted sign, to use in exp(...)
        final double total_dE =   // Calculate potential change in fromBox
                potentialOfPair(fromBox,
                                trialIndex,
                                new double[]{prtcls[0][0][trialIndex], prtcls[1][0][trialIndex]},
                                new double[]{prtcls[0][1][trialIndex], prtcls[1][1][trialIndex]},
                                new double[]{prtcls[0][2][trialIndex], prtcls[1][2][trialIndex]})
                        -
                        potentialOfPair(toBox,
                                        trialIndex,
                                        xTrial,
                                        yTrial,
                                        zTrial); // same change in toBox


        // XXX todo: remove after debug
//        try {
//            System.out.println("delta E: " + total_dE);
//            Thread.sleep(700);
//        } catch (InterruptedException e1) {
//            e1.printStackTrace();
//        }

        return FastMath.exp(total_dE) * expoCoefficient;
    }

    /**
     * calculates the potential energy value of a given pair inside a box (in the main 'prtcls'
     * array)
     *
     * @param box       box ralative to which potential is calculated
     * @param pairIndex which pair inside prtcls to calculate
     * @param x         double[] {[type1], [type2]} – x coordinates of the pair (either current or test)
     * @param y         double[] {[type1], [type2]} – y coordinates of the pair (either current or test)
     * @param z         double[] {[type1], [type2]} – z coordinates of the pair (either current or test)
     * @return the value of a potential energy that given pair would add up to the total box energy
     */
    private final double potentialOfPair(int box, int pairIndex,
                                         final double[] x,
                                         final double[] y,
                                         final double[] z
    ) {
        return sumPotential(0, pairIndex, // first particle energy in old box
                            x[0], y[0], z[0],
                            prtcls,
                            box, halfBox[box]) +
                sumPotential(1, pairIndex, // second particle energy in old box
                             x[1], y[1], z[1],
                             prtcls,
                             box, halfBox[box])
                - // minus the interparticle energy (counted twice)
                getPotential(0,
                             x[0], y[0], z[0],
                             1,
                             x[1], y[1], z[1],
                             halfBox[box]);
    }

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

        xTrial[0] = correctPosition(prtcls[trialType][0][trialIndex] + x, boxSize[lastBox]);
        yTrial[0] = correctPosition(prtcls[trialType][1][trialIndex] + y, boxSize[lastBox]);
        zTrial[0] = correctPosition(prtcls[trialType][2][trialIndex] + z, boxSize[lastBox]);

        // Calculating the potential shift
        // firstly, old Energy
        final double oldE = sumPotential(trialType, trialIndex,
                                         prtcls[trialType][0][trialIndex],
                                         prtcls[trialType][1][trialIndex],
                                         prtcls[trialType][2][trialIndex],
                                         prtcls,
                                         lastBox, halfBox[lastBox]
        );

        // then, new Energy
        final double newE = sumPotential(trialType,
                                         trialIndex,
                                         xTrial[0],
                                         yTrial[0],
                                         zTrial[0],
                                         prtcls,
                                         lastBox, halfBox[lastBox]
        );


        return newE - oldE;
    }

    /**
     * Sum of the potentials between given particle (with given coords) and all other particles in
     * the box
     *
     * @param whichBox in which box is summation
     * @return
     */
    private double sumPotential(
            final int particleType, final int particleIndex,
            final double x, final double y, final double z,
            final double[][][] particls,
            final int whichBox,
            final double halfbox
    ) {
        final int length = lengths[whichBox];
        final int offset = whichBox * boxBorder;

        double result = 0.0;


        for (int type = 0; type < 2; type++) {
            for (int i = offset; i < length + offset; i++) {
                if (type != particleType || i != particleIndex) {
                    result = result + getPotential(type,
                                                   particls[type][0][i],
                                                   particls[type][1][i],
                                                   particls[type][2][i],
                                                   particleType,
                                                   x, y, z, halfbox);
                }
            }
        }
        return result;
    }


    private final double dSquared(double dx, double dy, double dz, final double halfbox) {
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

    /*
     * Two particle potential value between given particles with given coordinates.
     *
     * @param type1,2  – first particle type in prtcls[type][][] array
     * @param num1     – first particle number in prtcls[][][num] array
     * @param particls - which particles array to use (either main or test)
     * @param whichBox - which box parameters to use (boxSize, halfSize, etc.)
     * @return energy value in kT units
     */
    /*private final double getPotential(
            final int type1, final int num1,
            final int type2, final int num2,
            final double trialX, final double trialY, final double trialZ,
            final double[][][] particls,
            final int whichBox
    ) {
        if (type1 == type2 && num1 == num2) {
            return 0;
        }

        return getPotential(type1,
                            particls[type1][0][num1],
                            particls[type1][1][num1],
                            particls[type1][2][num1],
                            type2, trialX, trialY, trialZ,
                            whichBox);
    } */

    private final double getPotential(final int type1,
                                      final double x1,
                                      final double y1,
                                      final double z1,
                                      final int type2,
                                      final double x2,
                                      final double y2,
                                      final double z2,
                                      final double halfbox
    ) {
        final double r = Math.sqrt(dSquared(x1 - x2, y1 - y2, z1 - z2, halfbox));

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
    @Override
    public double[] getCurrentResult() {
        return new double[]{
                gamma(densitiesAvg[0]),
                reducedEnrgyAvg[0],
                gamma(densitiesAvg[1]),
                reducedEnrgyAvg[1],
        };
    }

    /**
     * @return gamma for given density in cm^-3. For current ensemble, using its T value
     */
    private final double gamma(double density) {return e * e * FastMath.cbrt(density) / (k * T);}

    /**
     * @param box for which box
     * @return current reduced volume for given box (1/gamma^3) in the ensemble
     */
    protected final double getVstar(int box) {
        final double onegamma = (k * T) / (e * e);
        return onegamma * onegamma * onegamma / densitiesAvg[box];
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

    void setCurrEnergies(double en1, double en2) {
        currentEnergy[0] = en1;
        currentEnergy[1] = en2;
    }

    double[] getDensitiesAvg() {
        return densitiesAvg;
    }

    void setDensitiesAvg(final double dens0, final double dens1) {
        densitiesAvg[0] = dens0;
        densitiesAvg[1] = dens1;
    }

    double[][][] getParticles() {return prtcls;}

    public double getAvgEnergy(int box) {
        return reducedEnrgyAvg[box];
    }
}
