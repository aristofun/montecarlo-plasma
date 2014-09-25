package com.butlitsky.mk.ensembles;

import com.butlitsky.mk.IEnsemble;
import com.butlitsky.mk.math.MersenneTwisterFast;
import com.butlitsky.mk.options.CLOptions;
import com.butlitsky.mk.options.EOptions;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.concurrent.ThreadLocalRandom;

/**
 * Metropolis algorithm implementation suitable for any ensemble type.
 * <p/>
 * Created by aristofun on 21.07.14.
 */
public abstract class MetropolisEnsemble implements IEnsemble {

    /* Just number formatting utilities */
    protected final NumberFormat FORMAT = new DecimalFormat(EOptions.SCIENTIFIC_FORMAT_STR);
    protected final NumberFormat SHORT_FORMAT = new DecimalFormat(EOptions.SHORT_FORMAT_STR);
    protected final NumberFormat MICRO_FORMAT = new DecimalFormat(EOptions.MICRO_FORMAT_STR);
    /**
     * Every Nth step to calculate rare or heavy values (e.g. save configuration, correlation etc.)
     * ~ numPart * 7
     */
    private final int CALC_RARE_INT;

    /**
     * Every Nth step to calculate occasionally needed values (e.g. correlation array etc.)
     * ~ numPart * 3
     */
    private final int CALC_MID_INT;

    /**
     * Every Nth step to calculate frequently needed values (e.g. energy averages)
     * ~ numPart by default
     */
    private final int CALC_FREQUENT_INT;

    /**
     * Options object for the ensemble
     */
    protected final EOptions opt;

    /**
     * Folder to store/load results and states to/from.
     */
    protected final String myFolder;

    /**
     * Should be overriden by subclass to indicate the ensemble point if need more details.
     * Default is "7K/2e12" format
     */
    protected String myTag;

    /**
     * running control flag
     */
    private boolean finished = false;


    // ------------ Monte Karlo --------------
    /**
     * total number of steps to run
     */
    private final int numSteps;

    /**
     * currently running step
     */
    private int currStep;

    /**
     * Hi quality random generator
     */
    private final MersenneTwisterFast rnd;

    /**
     * simple random generator (for choosing particle for example, etc.)
     */
    private final ThreadLocalRandom localRnd;


    protected MetropolisEnsemble(EOptions options, int frequentInterval, int midInterval,
                                 int rareInterval) {
        rnd = new MersenneTwisterFast();
        localRnd = ThreadLocalRandom.current();

        opt = options;
        myFolder = options.getFolder();
        numSteps = options.getNumSteps();
        myTag = options.getT() + "/" + MICRO_FORMAT.format(options.getDensity());

        CALC_FREQUENT_INT = frequentInterval;
        CALC_MID_INT = midInterval;
        CALC_RARE_INT = rareInterval;
    }


    @Override
    public int getCurrStep() {
        return currStep;
    }

    @Override
    public int getNumSteps() {
        return numSteps;
    }

    protected int nextInt(int numPart) {
        return localRnd.nextInt(numPart);
    }

    @Override
    public boolean isFinished() {
        return finished;
    }

    @Override
    public void stop() {
        finished = true;
    }

    /**
     * All subclasses must be sure to load their states before running the execution. I.e.
     * loadState() must be called before running.
     */
    @Override
    public void run() {
        if (currStep >= numSteps || finished) {
//            System.out.print(myFolder + " No run.\t");
            finished = true;
            return;
        }

        int i = currStep;

        // fast run through initial steps, unstable configuration
        if (currStep < CLOptions.INITIAL_STEPS && numSteps > CLOptions.INITIAL_STEPS) {
            System.out.println("Ignoring first " + CLOptions.INITIAL_STEPS + " steps...");

            while (i < CLOptions.INITIAL_STEPS) {
                if (finished) {
                    System.out.println("STOP " + myFolder + ", finished=true\t");
                    break;
                }
                play(i);
                i++;
            }
        }

        currStep = i;

        if (!finished) {
            while (i < numSteps) {
                if (finished) {
                    System.out.println("STOP " + myFolder + ", finished=true\t");
                    break;
                }

                if (play(i)) {
                    onTrialAccepted();
                    currStep = i;
                } else {
                    onTrialRejected();
                }

                if (i % CALC_FREQUENT_INT == 0) {
                    doFrequentCalc();
                }
                if (i % CALC_MID_INT == 0) {
                    doMidCalc();
                }
                if (i % CALC_RARE_INT == 0) {
                    doRareCalc();
                }

                i++;
            }

        }
        currStep = i;
        finished = true;

        saveStateOnStop();

        System.out.print("" + myFolder + " finished on " + currStep + " steps.\t");
    }


    /**
     * Main Metropolis method.
     * It should make a trial random change and check the criteria if it should be accepted or not.
     * <p/>
     * This method MUST change the current ensemble real state (if the trial move is accepted)!
     * <p/>
     * So that calling this method sequentially eventually changes the ensemble state to more
     * probable.
     *
     * @param step
     * @return true &mdash; if the trial move was accepted, false &mdash; otherwise.
     */
    protected abstract boolean play(int step);


    /**
     * Do what you have to do after trial move is rejected.
     */
    protected abstract void onTrialRejected();

    /**
     * Do what you have to do after trial move is accepted.
     */
    protected abstract void onTrialAccepted();

    /**
     * Every ~numPart*7 steps action
     * <p/>
     * // saveCorrelation();
     * // saveState();
     */
    protected abstract void doRareCalc();

    /**
     * Every ~numPart*3 steps action
     * <p/>
     * //   averageEnergy();
     * // if (saveLongTail) saveLongTail();
     * //calcCorrelation();
     */
    protected abstract void doMidCalc();

    /**
     * Every ~numPart steps action
     */
    protected abstract void doFrequentCalc();

    /**
     * random (-1.0; 1.0)
     */
    final double myRandom() {
        return rnd.nextDouble(true, true) * 2.0 - 1.0;
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
     * Contract method to save results of a Metropolis calculation.
     * Called at least in the end of Metropolis chain or if the Ensemble is gracefully stopped
     * using {@link #stop()} method.
     * <p/>
     * It's up to children classes to call it during the calculation to save intermediate state.
     * <p/>
     * //        saveLongTail();
     * //        closeLongTail();
     *
     * @see #doFrequentCalc()
     * @see #doMidCalc()
     * @see #doRareCalc()
     */
    protected abstract void saveStateOnStop();

    @Override
    public abstract void loadState();

    /**
     * Used to set current step to run from. After recovering state from a file.
     *
     * @param step
     */
    protected void setCurrStep(int step) {
        currStep = step;
    }


    /**
     * Short summary to distinct one ensemble from another in the interface.
     * <p/>
     * Warning! Must be use for info purposes only! For correct files & states manipulation use
     * {@link com.butlitsky.mk.options.EOptions#getFolder()} of the Ensemble.
     *
     * @return string tag, for example "1K/2e7" for current ensemble point
     */
    @Override
    public String toString() {
        return myTag;
    }
}
