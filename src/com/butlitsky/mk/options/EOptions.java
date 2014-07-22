package com.butlitsky.mk.options;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.LinkedHashSet;
import java.util.Scanner;
import java.util.Set;

import static com.butlitsky.mk.IEnsemble.e;
import static com.butlitsky.mk.IEnsemble.k;
import static java.lang.Math.abs;
import static java.lang.Math.pow;

/**
 * Specific options instance for each calculation point. Initialized by config file + CLOptions
 */
public class EOptions implements CharSequence {
    public static final String SCIENTIFIC_FORMAT_STR = "0.000000000000000E0";
    public static final String SHORT_FORMAT_STR = "0.000E0";
    public static final String MICRO_FORMAT_STR = "0.##E0";

    private final NumberFormat MICRO_FORMAT = new DecimalFormat(MICRO_FORMAT_STR);


    private final double myDensity;
    private final double maxDelta;
    private final int myNumParticles;
    private final int myNumSteps;

    /**
     * STRATEGY bits 0 – default, 1 – save longtail
     */
    private final int myStrategy;

    private final int T;

    private boolean isOld;

    private final double gamma;

    private final String myFolder;

    EOptions(int t, double density, double dX, int numParticles,
             int numSteps, int strat, boolean old) {
        myDensity = density;
        maxDelta = dX;
        myNumParticles = numParticles;
        myNumSteps = numSteps;
        myStrategy = strat;
        T = t;
        isOld = old;

        gamma = e * e * pow(2 * density, 0.333333333333333) / (k * T);

        myFolder = "_" + T + "K_" + numParticles + "pa_d" + maxDelta + "/" +
                MICRO_FORMAT.format(myDensity);
    }

    public double getDensity() {
        return myDensity;
    }

    public double getMaxDelta() {
        return maxDelta;
    }

    public final int getNumParticles() {
        return myNumParticles;
    }

    public int getNumSteps() {
        return myNumSteps;
    }

    public int getT() {
        return T;
    }

    public int getStrategy() {
        return myStrategy;
    }


    public boolean isOld() {
        return isOld;
    }

    public void setOld(boolean flag) {
        isOld = flag;
    }

    /**
     * Two EOptions are equal if they point to the same folder where calculations states are stored
     */
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
        return myFolder;
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
        return T + " " + MICRO_FORMAT.format(myDensity) + " " + MICRO_FORMAT.format(maxDelta) + " " +
                myNumParticles + " " + myNumSteps + " " + myStrategy + " " + isOld;
    }

    public double getGamma() {
        return gamma;
    }


    /**
     * Create EOptions instance for current calculation point parameters. Some properties are
     * taken from {@link CLOptions}
     *
     * @param line Options format is: <code>T (in K), Density (Ne=Ni 1/cm3), numSteps,
     *             strategy (longtail saving - 1), isOld</code>
     */
    private static EOptions createFromLine(String line) {

        Scanner s = new Scanner(line);


        int t = s.nextInt();
        double density = s.nextDouble();
        int numSteps = s.nextInt();
        int strategy = s.nextInt();
        boolean isOld = s.nextBoolean();

        return new EOptions(
                (CLOptions.DEFAULT_TEMP < 0) ? abs(t) : CLOptions.DEFAULT_TEMP, // T
                abs(density),    // density
                CLOptions.MAX_DELTA_FACTOR, // delta
                CLOptions.NUM_PARTICLES, // numPart
                (CLOptions.DEFAULT_NUM_STEPS < 0) ? numSteps : CLOptions.DEFAULT_NUM_STEPS, // numSteps
                strategy,       // strategy
                isOld);       // isOld
    }

    /**
     * Reading configuration from file assuming line by line options configuration
     *
     * @param configFile
     * @return
     */
    public static Set<EOptions> readConfig(String configFile) throws IOException {
        Set<EOptions> opts = new LinkedHashSet<>();

        try {
            BufferedReader bufRead = new BufferedReader(new FileReader(configFile));

            String line;    // String that holds current file line
            // Read through file one line at time.
            do {
                line = bufRead.readLine();
                if (line != null && !line.trim().startsWith("#") && !line.trim().equals(""))
                    opts.add(createFromLine(line.trim()));

            } while (line != null);

            bufRead.close();

        } catch (FileNotFoundException e) {
            System.out.println("File " + configFile + " not found");
            throw e;
        } catch (IOException e) {
            System.out.println("Config file broken, i'm quit.");
            throw e;
        }
        opts.remove(null); // just in case of blank EnsembleOptions were added

        return opts;
    }
}
