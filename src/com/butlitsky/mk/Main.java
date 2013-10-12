package com.butlitsky.mk;

import org.apache.commons.cli.*;

import java.util.Date;
import java.util.Locale;

@SuppressWarnings("AccessStaticViaInstance")
public class Main {
    public static final String version = "6.0 _ asymmetric potentials & energy introduced";

    /**
     * <pre>usage: (./runmk.command | runmk.bat) [OPTIONS]
     * -d,--delta <DELTA_FACTOR>   maxDx coeff. (1.5 default)
     * -h                          show this help and exit
     * -pa,--particles <NUM>       number of particles, if set all mk_config.ini options ignored
     * -stp,--steps <NUM>          number of total steps, if set all mk_config.ini options ignored
     * -po,--polka <POLKA>         polochka parameter value (2.0 default)
     * -r,--refresh <SECONDS>      threads status refresh interval (30 sec. default)
     * -w,--workers <NUM>          number of parallel threads (default is MAX(2, CPUs/2)
     * </pre>
     */
    public static void main(String[] args) {

        Date start = new Date();
        System.out.println("\nMonte-Karlo game v. " + version + ", (c) Michael Butlitsky 2013\n");

        /* playground
        Deque<Double> deck = new ArrayDeque<Double>(3);
        deck.addFirst(4.9);
        deck.addFirst(3.9);
        deck.addFirst(0.9);

        while (start != null) {
            deck.addFirst(7.7);
            deck.addFirst(2.7);
            deck.removeLast();
            deck.removeLast();

            System.out.println(deck);
        }             */

        // apache CLI lib options parser
        parseArgs(args);

        System.out.println();
        System.out.println(start);// System.currentTimeMillis());
        Locale.setDefault(Locale.US); // for reading/writing '.' delimited doubles properly

        try {
            final EnsembleController controller = new EnsembleController();

            Runtime.getRuntime().addShutdownHook(new Thread() {
                public void run() {
                    controller.stop();
                    try {
                        Thread.sleep(EnsembleController.REFRESH_DELAY);
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    System.out.println("Got SHUTDOWN hook, gracefully complete.");
                }
            });


//            controller.saveContinueOptions();
            controller.start();


        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("FATAL: something fucked up – " + e.getLocalizedMessage());
            System.exit(1);
        }

        Date fin = new Date();
        System.out.println(fin);
        System.out.println("\nJob took " + getHumanTimeDiff(start, fin) + ". Bye!\n");
    }

    private static void parseArgs(String[] args) {
//        Option temp = OptionBuilder.withArgName("TEMP").hasArg().withDescription("temperature " +
//                "for polochka (1K default)").withLongOpt("temp").create("t");

        Option polka = OptionBuilder.withArgName("POLKA").hasArg().withDescription("polochka " +
                "parameter value (2.0 default)").withLongOpt("polka").create("po");

        Option workers = OptionBuilder.withArgName("NUM").hasArg().withDescription("number of " +
                "parallel threads (default is MAX(2, CPUs/2)").withLongOpt("workers").create("w");

        Option particles = OptionBuilder.withArgName("NUM").hasArg().withDescription("number of " +
                "particles, if set all mk_config.ini options ignored").withLongOpt("particles")
                .create("pa");

        Option avpoints = OptionBuilder.withArgName("NUM").hasArg().withDescription("number of " +
                "averaging points for Energy (128 default)").withLongOpt("avpoints")
                .create("ap");

        Option steps = OptionBuilder.withArgName("NUM").hasArg().withDescription("number of " +
                "steps (x Number of Particles), if set all mk_config.ini options ignored")
                .withLongOpt("steps")
                .create("stp");

        Option delta = OptionBuilder.withArgName("DELTA_FACTOR").hasArg()
                .withDescription("maxDx coeff. (overrides .ini parameters if set) \n" +
                        "        – zero equals 1. x BOX SIZE \n" +
                        "        – if float number >= 1  trial shift delta position == this " +
                        "number * average distance between particles (depends on density) along every axis.\n" +
                        "        – if float number < 1 trial shift delta is this fraction of a " +
                        "box size along each axis (i. e. 0.5 means trial particle shifts to half of box size along each axis). \n")
                .withLongOpt("delta").create("d");

        Option refresh = OptionBuilder.withArgName("SECONDS").hasArg()
                .withDescription("threads status refresh interval (5 sec. default)")
                .withLongOpt("refresh").create("r");

        Option ewaldDelta = OptionBuilder.withArgName("NUM").hasArg().withDescription("Ewald " +
                "accuracy delta parameter (0.001 default)").withLongOpt("ewaldelta").create("ewd");

        Option ewaldNcut = OptionBuilder.withArgName("NUM").hasArg().withDescription("Ewald" +
                "Ncutoff parameter (3 default)").withLongOpt("ewaldn").create("ewn");

        Option harrisR = OptionBuilder.withArgName("NUM").hasArg().withDescription("Harris " +
                "radius. When set Harris algorithm used with given R in L (5 default)")
                .withLongOpt("harris").create("harris");

        Options options = new Options();
        options.addOption("h", false, "show this help and exit");
        options.addOption("ew", false, "use Ewald summation");
        options.addOption(avpoints);
        options.addOption(polka);
        options.addOption(delta);
        options.addOption(refresh);
        options.addOption(workers);
        options.addOption(particles);
        options.addOption(steps);
        options.addOption(ewaldDelta);
        options.addOption(ewaldNcut);
        options.addOption(harrisR);

        CommandLineParser parser = new BasicParser();

        try {
            // parse the command line arguments
            CommandLine line = parser.parse(options, args);

            if (line.hasOption("h")) {
                // automatically generate the help statement
                printHelpAndExit(options);
            }

            if (line.hasOption("ew")) {
                EnsembleController.ENS_TYPE = 1;
                System.out.println("EWALD calculations!");

                if (line.hasOption("ewn")) {
                    EnsemblePolochkaEwald.Ncutoff = Integer.parseInt(line.getOptionValue("ewn"));
                }
                if (line.hasOption("ewd")) {
                    EnsemblePolochkaEwald.DELTA = Double.parseDouble(line.getOptionValue("ewd"));
                }
            }

            if (line.hasOption("harris")) {
                EnsembleController.ENS_TYPE = 2;
                System.out.println("Harrison calculations!")
                ;
                EnsemblePolochkaHarrison.N
                        = Integer.parseInt(line.getOptionValue("harris"));
                System.out.println("Harrison R = " + EnsemblePolochkaHarrison.N + "");
            }


            if (line.hasOption("pa")) {
                Ensemble.DEFAULT_NUM_PARTICLES = Integer.parseInt(line.getOptionValue("pa"));
                System.out.println("Custom particles number = " + Ensemble.DEFAULT_NUM_PARTICLES);
            }

            if (line.hasOption("stp")) {
                Ensemble.DEFAULT_NUM_STEPS = Integer.parseInt(line.getOptionValue("stp"));
                System.out.println("Custom steps number = " + Ensemble.DEFAULT_NUM_STEPS);
            }

            if (line.hasOption("po")) {
                EnsemblePolochka.EPSILON = Double.parseDouble(line.getOptionValue("polka"));
            }
            System.out.println("Polochka (electron-ion) = - " + EnsemblePolochka.EPSILON);

//            if (line.hasOption("t")) {
//                EnsemblePolochka.EPSILON = Double.parseDouble(line.getOptionValue("polka"));
//            }
//            System.out.println("Polochka (electron-ion) = – " + EnsemblePolochka.EPSILON);

            if (line.hasOption("w")) {
                EnsembleController.NUM_THREADS = Integer.parseInt(line.getOptionValue("workers"));
                System.out.println("Custom workers count = " + EnsembleController.NUM_THREADS);
            }

            if (line.hasOption("d")) {
                Ensemble.DELTA_FACTOR = Double.parseDouble(line.getOptionValue("delta"));
                System.out.println("Custom delta X = " + Ensemble.DELTA_FACTOR);
            }

            if (line.hasOption("r")) {
                EnsembleController.REFRESH_DELAY
                        = 1000 * Integer.parseInt(line.getOptionValue("refresh"));
            }
            System.out.println("Status refresh delay = " + EnsembleController.REFRESH_DELAY / 1000);

            if (line.hasOption("ap")) {
                Ensemble.NUM_ENERGY_AVG_POINTS
                        = Integer.parseInt(line.getOptionValue("ap"));
            }
            System.out.println("Avg. points = " + Ensemble.NUM_ENERGY_AVG_POINTS);
        } catch (Exception exp) {
            // oops, something went wrong
            System.err.println("Parsing failed, using defaults  Reason: " + exp.getMessage());
            printHelpAndExit(options);
        }
    }

    private static void printHelpAndExit(Options options) {
        HelpFormatter formatter = new HelpFormatter();
        formatter.printHelp("(./runmk.command | runmk.bat) [OPTIONS]", options);
        System.exit(0);
    }

    private static String getHumanTimeDiff(Date oldTime, Date newTime) {
        long diffInSeconds = (newTime.getTime() - oldTime.getTime()) / 1000;

        long diff[] = new long[]{0, 0, 0, 0};
    /* sec */
        diff[3] = (diffInSeconds >= 60 ? diffInSeconds % 60 : diffInSeconds);
    /* min */
        diff[2] = (diffInSeconds = (diffInSeconds / 60)) >= 60 ? diffInSeconds % 60 : diffInSeconds;
    /* hours */
        diff[1] = (diffInSeconds = (diffInSeconds / 60)); // >= 24 ? diffInSeconds % 24 : diffInSeconds;
//    /* days */
//        diff[0] = (diffInSeconds = (diffInSeconds / 24));

        return String.format(
//                "%d day%s, %d hour%s, %d minute%s, %d second%s",
                "%d hour%s, %d minute%s, %d second%s",
//                diff[0],
//                diff[0] > 1 ? "s" : "",
                diff[1],
                diff[1] > 1 ? "s" : "",
                diff[2],
                diff[2] > 1 ? "s" : "",
                diff[3],
                diff[3] > 1 ? "s" : "");
    }
}
