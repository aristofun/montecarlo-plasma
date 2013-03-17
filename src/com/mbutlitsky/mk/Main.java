package com.mbutlitsky.mk;

import org.apache.commons.cli.*;

import java.util.Date;
import java.util.Locale;

@SuppressWarnings("AccessStaticViaInstance")
public class Main {
    public static final String version = "0.3";

    /**
     * @param args First – polochka deepness, Second – maxDelta factor times average distance
     *             defaults are 2.0 and 1.3
     */
    public static void main(String[] args) {
        Date start = new Date();
        System.out.println("Monte-Karlo game v. " + version + ", (c) Michael Butlitsky 2013");

        // apache CLI lib options parser
        parseArgs(args);

        System.out.println("Polochka set = " + EnsemblePolochka.EPSILON);
        System.out.println("Delta factor set = " + Ensemble.DELTA_FACTOR);

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


            controller.saveContinueOptions();
            controller.start();


        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("FATAL: something fucked up: " + e.getLocalizedMessage());
            System.exit(1);
        }

        Date fin = new Date();
        System.out.println(fin);
        System.out.println("Job took " + getHumanTimeDiff(start, fin) + ". Bye!");
    }

    private static void parseArgs(String[] args) {
        Option polka = OptionBuilder.withArgName("POLKA").hasArg().withDescription("polochka " +
                "parameter value (2.0 default)").create("polka");

        Option delta = OptionBuilder.withArgName("DELTA_FACTOR").hasArg()
                .withDescription("maxDx coeff. (1.3 default)").create("delta");

        Option refresh = OptionBuilder.withArgName("SECONDS").hasArg()
                .withDescription("threads status refresh interval (20 sec. default)")
                .create("refresh");


        Options options = new Options();
        options.addOption("h", false, "show this help and exit");
        options.addOption(polka);
        options.addOption(delta);
        options.addOption(refresh);

        CommandLineParser parser = new BasicParser();

        try {
            // parse the command line arguments
            CommandLine line = parser.parse(options, args);

            if (line.hasOption("h")) {
                // automatically generate the help statement
                printHelpAndExit(options);
            }

            if (line.hasOption("polka")) {
                EnsemblePolochka.EPSILON = Double.parseDouble(line.getOptionValue("polka"));
            }

            if (line.hasOption("delta")) {
                Ensemble.DELTA_FACTOR = Double.parseDouble(line.getOptionValue("delta"));
            }

            if (line.hasOption("refresh")) {
                EnsembleController.REFRESH_DELAY
                        = 1000 * Integer.parseInt(line.getOptionValue("refresh"));
            }

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
