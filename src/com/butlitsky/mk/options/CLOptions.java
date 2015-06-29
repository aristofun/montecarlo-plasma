package com.butlitsky.mk.options;

import org.apache.commons.cli.*;

/**
 * Global options storage with CL arguments parser.
 * <p/>
 * It knows nothing about implementation except the options it stores.
 * <p/>
 * Created by aristofun on 21.07.14.
 */
public class CLOptions {

    /**
     * 0 – default Polochka, no long range account
     * 1 – basic Ewald summation
     * 2 – Harrison algorithm
     * 3 – Pseudo Potential
     * 4 – Gibbs ensemble (two separate boxes)
     * 5 – Gibbs ensemble (two separate boxes) Lennard-Johnes
     */
    public static int ENSEMBLE_TYPE = 0;

    /**
     * how often check child ensembles & display status on console (in msec)
     */
    public static int REFRESH_DELAY = 15000;

    /**
     * may be overriden by CLI options, if zero – max(2, CPUs/2) used
     */
    public static int NUM_THREADS = -1;


// ----------------- Common ensemble options --------------------------------
    /**
     * overrides specific particle number in ini file if set in CLI options
     */
    public static int NUM_PARTICLES = 500;

    /**
     * number of steps for all ensembles
     */
    public static int DEFAULT_NUM_STEPS = -1;

    /**
     * how many steps from beginning to ignore  in Markov chain
     */
    public static int INITIAL_STEPS = 50000;

    /**
     * Default polochka deepness in kT (Epsilon parameter for Lennard-Johnes case = 1/Tstar)
     */
    public static double POLOCHKA = 4.0;

    /**
     * temperature to be overriden by command line or *.ini options
     */
    public static int DEFAULT_TEMP = -1;

    /**
     * These many recent Metropolis steps is taken into account for averages calculation
     */
    public static int NUM_ENERGY_AVG_STEPS = 100000;

    /**
     * Maximum dX,dY,dZ random displacement coefficient (times avg. distance)
     */
    public static double MAX_DELTA_X = 1;

    public static boolean START_FROM_FCC = false;

// -------------  Ewald ensemble specific steps --------------------------
    /**
     * Ewal n cutoff
     */
    public static int EWALD_N_CUTOFF = 3;

    /**
     * Ewald alpha
     */
    public static double EWALD_DELTA = 0.00000001;

// -------------  Harrison ensemble specific steps -----------------------
    /**
     * R = N * BoxSize
     */
    public static int HARRISON_N = 1;

// -------------- Gibbs ensemble specific options ------------------------
    /**
     * current ensembles values resolution by steps (used to print data for plotting
     * the dynamics of ensembles states migration).
     */
    public static int N_RESOLUTION_STEPS = 1000;

    /**
     * Initial Ro* parameter for LJ ensemble only
     */
    public static double RO_STAR = 0.1;

    /**
     * Gibbs maximum relative volume change from 0 to 1
     */
    public static double MAX_DELTA_V = 0.07;


    private CLOptions() {
    }

    public static String getOneLineSummary() {
        return
                "Ens.type " + ENSEMBLE_TYPE + ", Num.part " + NUM_PARTICLES + ", " +
                        "steps " + DEFAULT_NUM_STEPS / 1000
                        + "K, polka " + POLOCHKA + ", dX " + MAX_DELTA_X + ", Avg.points "
                        + NUM_ENERGY_AVG_STEPS / 1000 + "K, refresh " + REFRESH_DELAY / 1000;
    }

    public static void init(String[] args) throws IllegalArgumentException {
        Options options = buildOptions();

        CommandLineParser parser = new BasicParser();
        CommandLine line = null;

        try {
            line = parser.parse(options, args);
        } catch (ParseException e) {
            printHelpAndExit(options);
            throw new IllegalArgumentException();
        }


        if (line.hasOption("h")) {
            // automatically generate the help statement
            printHelpAndExit(options);
            throw new IllegalArgumentException();
        }

        if (line.hasOption("cubic")) { // new in 8.0
            START_FROM_FCC = true;
            System.out.println("Initial NaCl fcc cubic");
        }

        if (line.hasOption("pseudo")) {
            ENSEMBLE_TYPE = 3;
            System.out.println("PSEUDO potential");
        }

        if (line.hasOption("gibbs")) { // new in v. 9.0
            ENSEMBLE_TYPE = 4;
            System.out.println("Gibbs ensemble");
        }

        if (line.hasOption("gibbs_lj")) { // new in v. 11.0
            ENSEMBLE_TYPE = 5;
            System.out.println("Gibbs LJ ensemble");

            if (line.hasOption("rostar")) {
                RO_STAR = Double.parseDouble(line.getOptionValue("rostar"));
            }
        }

        if (line.hasOption("ew")) {
            ENSEMBLE_TYPE = 1;
            System.out.println("EWALD calculations!");

            if (line.hasOption("ewn")) {
                EWALD_N_CUTOFF = Integer.parseInt(line.getOptionValue("ewn"));
            }
            if (line.hasOption("ewd")) {
                EWALD_DELTA = Double.parseDouble(line.getOptionValue("ewd"));
            }
        }

        if (line.hasOption("harris")) {
            ENSEMBLE_TYPE = 2;
            System.out.println("Harrison calculations!");

            HARRISON_N = Integer.parseInt(line.getOptionValue("harris"));
            System.out.println("Harrison R = " + HARRISON_N + "");
        }


        if (line.hasOption("pa")) {
            NUM_PARTICLES = Integer.parseInt(line.getOptionValue("pa"));
        }

        System.out.println("Particles number = " + NUM_PARTICLES);

        if (line.hasOption("stp")) {
            DEFAULT_NUM_STEPS = Integer.parseInt(line.getOptionValue("stp"));
        }
        System.out.println("Steps number = " + DEFAULT_NUM_STEPS);

        if (line.hasOption("po")) {
            POLOCHKA = Double.parseDouble(line.getOptionValue("polka"));
        }
        System.out.println("Polochka (electron-ion) = - " + POLOCHKA);

        if (line.hasOption("t")) {
            DEFAULT_TEMP = Integer.parseInt(line.getOptionValue("temp"));
        }

        if (line.hasOption("w")) {
            NUM_THREADS = Integer.parseInt(line.getOptionValue("workers"));
            System.out.println("Custom workers count = " + NUM_THREADS);
        }

        if (line.hasOption("d")) {
            MAX_DELTA_X = Double.parseDouble(line.getOptionValue("delta"));
        }
        System.out.println("dX factor = " + MAX_DELTA_X);


        if (line.hasOption("inisteps")) { // new in 8.0
            INITIAL_STEPS = Integer.parseInt(line.getOptionValue("inisteps"));
        }
        System.out.println("INITIAL STEPS = " + INITIAL_STEPS);

        if (line.hasOption("r")) {
            REFRESH_DELAY = 1000 * Integer.parseInt(line.getOptionValue("refresh"));
        }
        System.out.println("Status refresh delay = " + REFRESH_DELAY / 1000);

        if (line.hasOption("ap")) {
            NUM_ENERGY_AVG_STEPS = Integer.parseInt(line.getOptionValue("ap"));
        }
        System.out.println("Avg. points = " + NUM_ENERGY_AVG_STEPS);

        if (line.hasOption("res")) {
            N_RESOLUTION_STEPS = Integer.parseInt(line.getOptionValue("res"));
        }
        System.out.println("N resolution = " + N_RESOLUTION_STEPS);
    }

    @SuppressWarnings("AccessStaticViaInstance")
    private static Options buildOptions() {
        Option temp = OptionBuilder.withArgName("TEMP").hasArg().withDescription("temperature " +
                                                                                         "for any ensemble type (overrides mk_config.ini, no default value)")
                .withLongOpt("temp").create("t");

        Option polka = OptionBuilder.withArgName("POLKA").hasArg().withDescription("polochka " +
                                                                                           "parameter value (" + POLOCHKA + " default)").withLongOpt(
                "polka").create("po");

        Option workers = OptionBuilder.withArgName("CPUs").hasArg().withDescription("number of " +
                                                                                            "parallel threads (default is MAX(2, CPUs/2)").withLongOpt(
                "workers").create("w");

        Option particles = OptionBuilder.withArgName("PARTICLES").hasArg().withDescription
                ("number of particles, (" + NUM_PARTICLES + " default)").withLongOpt("particles").create("pa");

        Option avpoints = OptionBuilder.withArgName("AVG.").hasArg().withDescription("number of " +
                                                                                             "averaging points for Energy (default is number of total workingssteps!)")
                .withLongOpt("avpoints").create("ap");


        Option steps = OptionBuilder.withArgName("STEPS").hasArg().withDescription("number of " +
                                                                                           "steps (" + DEFAULT_NUM_STEPS + " by default)").withLongOpt(
                "steps").create("stp");

        Option delta = OptionBuilder.withArgName("MAX_dX").hasArg()
                .withDescription("maxDx coeff. number * average distance between particles " +
                                         "(depends on density) along every axis (" + MAX_DELTA_X + " default)")
                .withLongOpt("delta").create("d");


        Option refresh = OptionBuilder.withArgName("SECONDS").hasArg()
                .withDescription("threads status refresh interval (" + REFRESH_DELAY / 1000 +
                                         "sec. default)").withLongOpt("refresh").create("r");

        Option ewaldDelta = OptionBuilder.withArgName("NUM").hasArg().withDescription("Ewald accuracy delta parameter ("
                                                                                              + EWALD_DELTA + " default)").withLongOpt(
                "ewaldelta").create("ewd");

        Option ewaldNcut = OptionBuilder.withArgName("NUM").hasArg().withDescription("Ewald cutoff parameter (" + EWALD_N_CUTOFF + " default)").withLongOpt(
                "ewaldn").create("ewn");

        Option harrisR = OptionBuilder.withArgName("NUM").hasArg().withDescription(
                "Harris radius. When set Harris algorithm used with given R in L (" + HARRISON_N + " default)")
                .withLongOpt("harris").create("harris");

        Option stepsToPass = OptionBuilder.withArgName("INI_STEPS").hasArg().withDescription
                ("Number of steps to ignore in markov chain averages (default " + INITIAL_STEPS + ")")
                .withLongOpt("inisteps").create("inisteps");

//        Gibbs specific options
        Option deltav = OptionBuilder.withArgName("MAX_dV/V").hasArg()
                .withDescription(" Gibbs maximum relative volume change from 0 to 1 (default is " +
                                         MAX_DELTA_V + " )").withLongOpt("deltav").create("dv");

        Option rostar = OptionBuilder.withArgName("RO_STAR").hasArg()
                .withDescription("Initial Lennard-Johnes Ro* parameter (0.1 default)")
                .withLongOpt("rostar").create("rostar");

        Option nResolution = OptionBuilder.withArgName("NUM").hasArg().withDescription(
                "number of averaging points for current values plotting (default is " + N_RESOLUTION_STEPS + " steps)")
                .withLongOpt("resolution").create("res");

        Options options = new Options();


//      common options
        options.addOption("h", false, "show this help and exit");
        options.addOption("pseudo", false, "use Pseudopotential");
        options.addOption("cubic", false, "use simple cubic start config (randomized by default)");
        options.addOption(avpoints);
        options.addOption(polka);
        options.addOption(delta);
        options.addOption(particles);
        options.addOption(refresh);
        options.addOption(workers);
        options.addOption(steps);
        options.addOption(temp);
        options.addOption(stepsToPass);

//      Gibbs options
        options.addOption("gibbs", false, "use Gibbs ensemble (two box simulation)");
        options.addOption("gibbs_lj", false, "Use gibbse ensemble calculation for Lennard-Johnes "); // new in v. 11.0
        options.addOption(deltav);
        options.addOption(nResolution);
        options.addOption(rostar);

//      Ewald options
        options.addOption("ew", false, "use Ewald summation");
        options.addOption(ewaldDelta);
        options.addOption(ewaldNcut);

//      Harrison options
        options.addOption(harrisR);

        return options;
    }

    private static void printHelpAndExit(Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("(./runmk.command | runmk.bat) [OPTIONS]", options);
    }

}
