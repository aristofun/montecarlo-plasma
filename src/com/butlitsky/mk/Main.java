package com.butlitsky.mk;

import com.butlitsky.mk.options.CLOptions;

import java.util.Date;
import java.util.Locale;

@SuppressWarnings("AccessStaticViaInstance")
public class Main {
    public static final String version = "13.0 / non-uniform LJ2 ensemble steps types";

    /**
     * <pre>
     * usage: (./runmk.command | runmk.bat) [OPTIONS]
     *  -ew                         use Ewald summation
     *  -ewd,--ewaldelta <NUM>      Ewald accuracy delta parameter (0.001 default)
     *  -ewn,--ewaldn <NUM>         EwaldNcutoff parameter (3 default)
     *
     *  -h                          show this help and exit
     *  -harris,--harris <NUM>      Harris radius. When set Harris algorithm used with given R in
     *  L (5 default)
     *
     *  -pseudo                     Use pseudo potential ensemble, all polochka parameters ignored
     *
     *  -ap,--avpoints <NUM>        number of averaging points for Energy (default is number of total steps!)
     *  -d,--delta <DELTA_FACTOR>   maxDx coeff.
     *      – zero equals  1. x BOX SIZE
     *      – other values = x avg. distance
     *
     *  -pa,--particles <NUM>       number of particles, if set all mk_config.ini options ignored
     *  -po,--polka <POLKA>         polochka parameter value (2.0 default)
     *  -
     *  -stp,--steps <NUM>          number of total steps to go for each point
     *  -r,--refresh <SECONDS>      threads status refresh interval (15 sec. default)
     *  -w,--workers <NUM>          number of parallel threads, default is MAX(2, CPUs/2)
     *  -cubic                      start from fcc NaCl initial configuration
     *
     *  –inisteps                   number of initial steps to ignore in averages calculations
     *
     *  –gibbs                      Use gibbse ensemble calculation (two boxes of total V and N)
     *  –gibbs_lj                   Use gibbse ensemble calculation for Lennard-Johnes potetial (two boxes of total V and N)
     *  –gibbs_lj2                   Use gibbse ensemble calculation for Lennard-Johnes potetial (two different boxes)
     *  –switch_rate                Percentage of interchange steps (default 0.05 - 5%)
     *  –rostar                     Initial Lennard-Johnes Ro* parameter (default 0.1)
     *  –rostar1                     Initial Lennard-Johnes Ro* parameter for first box (default 0.1)
     *  –rostar2                     Initial Lennard-Johnes Ro* parameter for second box (default 0.1)
     *  –N1                          Initial Lennard-Johnes N particles parameter for 1st box (default 0.1)
     *  –N2                          Initial Lennard-Johnes N parameter for 2nd box (default 0.1)
     *  –res, --resolution          Gibbs technique N step delta for per point current values plotting
     *  –dv, --deltav               Gibbs maximum relative volume change from 0 to 1 (default 0.15)
     * </pre>
     */
    public static void main(String[] args) {

        Date start = new Date();
        System.out.println("\nMonte-Karlo game v. " + version + ", (c) Michael Butlitsky 2013 + \n");

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
                        Thread.sleep(2000); // 2 sec to wait all finished
                    } catch (InterruptedException e) {
                        e.printStackTrace();
                    }
                    System.out.println("Got SHUTDOWN signal, gracefully complete.");
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
