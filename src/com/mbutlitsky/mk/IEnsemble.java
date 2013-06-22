package com.mbutlitsky.mk;

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
    public static final int CORR_LENGTH = 70;

    /**
     * current iteration step
     */
    public int getCurrStep();

    // Current state Parameters Getters & Setters ++++++++++

    double getEnergy();

//    double getPressure();

    int getNumPart();

    int getNumSteps();

    int getT();

    String getFolder();
    /**
     * returns true if this ensemble is not being executed anymore
     */
    boolean isFinished();

    /** gracefull stop the ensemble if running with saving state */
    void stop();
}

