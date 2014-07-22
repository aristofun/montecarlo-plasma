package com.butlitsky.mk;

/**
 * User: aristofun
 * Date: 02.03.13
 * Time: 16:58
 */
public interface IEnsemble extends Runnable {
    /**
     * in cm
     */
    public static final double BOHR = 5.2917721092e-9;
    /**
     * Elementary charge in SGS
     */
    public static final double e = 4.8032043e-10;

    /**
     * Boltsmann in SGS
     */
    public static final double k = 1.3806488e-16;

    public static final String STATE_FILE = "config.dat";
    public static final String LONGTAIL_FILE = "all_configs.dat";
    public static final String CORR_FILE = "correlation.dat";

    public static final int CORR_LENGTH = 90;
    /**
     * scaling factor to convert energy value to kT units
     * SCALE_FACTOR == e^2 / (Bohr * k)
     */
    double SCALE_FACTOR = 315775.01611746440408;

    /**
     * current iteration step
     */
    public int getCurrStep();

    // Current state Parameters Getters & Setters ++++++++++

    /**
     * Generic result for the ensemble in ready for use flavour (may be multiple values –
     * "energy, pressure" or "energy1, energy2" etc). Depends on the ensemble type implementation.
     */
    double[] getCurrentResult();

    /**
     * Contract method to load current state right before running Metropolis execution from
     * current step to numSteps.
     * <p/>
     * Must be called explicitly right after the instance creation (using factory pattern for
     * example).
     * <p/>
     * May either load state from a file or initialize.
     */
    void loadState();

    int getNumSteps();

    int getT();

    /**
     * returns true if this ensemble is not being executed anymore
     */
    boolean isFinished();

    /**
     * gracefull stop the ensemble if running with saving state
     */
    void stop();
}

