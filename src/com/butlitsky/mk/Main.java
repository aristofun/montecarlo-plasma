package com.butlitsky.mk;

import com.butlitsky.mk.options.CLOptions;

import java.util.Arrays;
import java.util.Date;
import java.util.Locale;

@SuppressWarnings("AccessStaticViaInstance")
public class Main {
    public static final String version = "9.0 + total refactoring of v.8";

    /**
     * <pre>
     * usage: (./runmk.command | runmk.bat) [OPTIONS]
     *  -ap,--avpoints <NUM>        number of averaging points for Energy (default is number of total steps!)
     *  -d,--delta <DELTA_FACTOR>   maxDx coeff. (overrides .ini parameters if set)
     *      – zero equals  1. x BOX SIZE
     *      – if float number >= 1 trial shift delta position == this number * average distance
     * between particles (depends on density) along every axis.
     *      – if float number < 1 trial shift delta is this fraction of a box size along each
     *        axis (i. e. 0.5 means trial particle shifts to half of box size along each axis).
     *  -ew                         use Ewald summation
     *  -ewd,--ewaldelta <NUM>      Ewald accuracy delta parameter (0.001 default)
     *  -ewn,--ewaldn <NUM>         EwaldNcutoff parameter (3 default)
     *  -h                          show this help and exit
     *  -harris,--harris <NUM>      Harris radius. When set Harris algorithm used with given R in
     *  L (5 default)
     *  -pseudo                     Use pseudo potential ensemble, all polochka parameters ignored
     *  -pa,--particles <NUM>       number of particles, if set all mk_config.ini options ignored
     *  -po,--polka <POLKA>         polochka parameter value (2.0 default)
     *  -r,--refresh <SECONDS>      threads status refresh interval (5 sec. default)
     *  -stp,--steps <NUM>          number of steps (x Number of Particles),
     *  if set all mk_config.ini options ignored
     *  -w,--workers <NUM>          number of parallel threads, default is MAX(2, CPUs/2)
     * </pre>
     */
    public static void main(String[] args) {

        Date start = new Date();
        System.out.println("\nMonte-Karlo game v. " + version + ", (c) Michael Butlitsky 2013\n");

        // apache CLI lib options parser


        parseArgs(args);

        System.out.println();
        System.out.println(start);// System.currentTimeMillis());
        Locale.setDefault(Locale.US); // for reading/writing '.' delimited doubles properly

        try {
            final IEnsembleController controller = new EnsembleController();

            Runtime.getRuntime().addShutdownHook(new Thread() {
                public void run() {
                    controller.stop();
                    try {
                        Thread.sleep(2000); // 1 sec to wait all finished
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    System.out.println("Got SHUTDOWN hook, gracefully complete.");
                }
            });

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
        try {
            // parse the command line arguments
            CLOptions.init(args);
        } catch (IllegalArgumentException exp) {
            // oops, not good options
            System.err.println("CL options not parsed.");
            System.exit(0);
        }
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
