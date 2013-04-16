package com.mbutlitsky.mk;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Scanner;

import static com.mbutlitsky.mk.IEnsemble.e;
import static com.mbutlitsky.mk.IEnsemble.k;
import static java.lang.Math.abs;
import static java.lang.Math.pow;

/**
 * Created with IntelliJ IDEA.
 * User: aristofun
 * Date: 02.03.13
 * Time: 17:10
 */
public class EOptions implements CharSequence {
    public static final String SCIENTIFIC_FORMAT_STR = "0.000000000000000E0";
    public static final String SHORT_FORMAT_STR = "0.0000E0";

    private final NumberFormat MICRO_FORMAT = new DecimalFormat("0.##E0");
    private final double density;
    private final double maxDelta;
    private final int numParticles;
    private final int numSteps;

    /** STRATEGY bits 0 – default, 1 – save longtail */
    private final int strategy;
    private final int T;
    private boolean isOld;
    private final double gamma;

    public EOptions(int t, double density, double maxDelta, int numParticles,
                    int numSteps, int strat, boolean old) {
        this.density = density;
        this.maxDelta = maxDelta;
        this.numParticles = numParticles;
        this.numSteps = numSteps;
        strategy = strat;
        T = t;
        isOld = old;

        gamma = e * e * pow(density, 0.333333333333333) / (k * T);
    }

    public static EOptions fromLine(String line) {
        Scanner s = new Scanner(line);
        //# T, Density, maxDx/y/z, numParticles, numSteps, strategy (not used yet), isOld
        return new EOptions(abs(s.nextInt()), abs(s.nextDouble()), abs(s.nextDouble()),
                abs(s.nextInt()), abs(s.nextInt()), abs(s.nextInt()), s.nextBoolean());
    }

    public double getDensity() {
        return density;
    }

    public double getMaxDelta() {
        return maxDelta;
    }

    public final int getNumParticles() {
        return numParticles;
    }

    public int getNumSteps() {
        return numSteps;
    }

    public int getT() {
        return T;
    }

    public int getStrategy() {
        return strategy;
    }


    public boolean isOld() {
        return isOld;
    }

    public void setOld(boolean flag) {
        isOld = flag;
    }

    /**
     * calculates gamma from N * T
     */
//    public double getGamma()  {
//
//    }
    public boolean equals(Object arg0) {
        if (this == arg0) {
            return true;
        }
        if (arg0 instanceof EOptions) {
            EOptions arg = (EOptions) arg0;
            return this.getFolder().equals(arg.getFolder());
        }
        return false;
    }

    public int hashCode() {
        return getFolder().hashCode();
    }

    /**
     * returns unique subfolder based on current options set
     */
    public String getFolder() {
        return T + "K/n_" + MICRO_FORMAT.format(density) + "";
    }

    @Override
    public int length() {
        return toString().length();
    }

    @Override
    public char charAt(int index) {
        return toString().charAt(index);
    }

    @Override
    public CharSequence subSequence(int start, int end) {
        return toString().subSequence(start, end);
    }

    public String toString() {
        return T + " " + MICRO_FORMAT.format(density) + " " + MICRO_FORMAT.format(maxDelta) + " " +
                numParticles + " " + numSteps + " " + strategy + " " + isOld;
    }

    public double getGamma() {
        return gamma;
    }

}
