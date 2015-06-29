package com.butlitsky.mk;

import com.butlitsky.mk.ensembles.EnsemblesFactory;
import com.butlitsky.mk.ensembles.GibbsConfigurationManager;
import com.butlitsky.mk.options.CLOptions;
import com.butlitsky.mk.options.EOptions;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadLocalRandom;

import static java.lang.Thread.sleep;

/**
 * Date: 02.03.13
 * Time: 19:07
 */
public class EnsembleController implements IEnsembleController {

    private volatile boolean running = true;

    // TODO: move this option to CLI
    private static final String CONFIG_FILE = "mk_config.ini";

    private final NumberFormat SHORT_FORMAT = new DecimalFormat(EOptions.SHORT_FORMAT_STR);
    private final NumberFormat FULL_FORMAT = new DecimalFormat(EOptions.SCIENTIFIC_FORMAT_STR);

    private final List<IEnsemble> ensembles = new ArrayList<>();

    private final Map<IEnsemble, Boolean> workerStates = new HashMap<>();

    private final Map<Integer, LinkedHashMap<IEnsemble, Deque<double[]>>> ensemblesResultValues =
            new HashMap<>();


    private final ExecutorService pool;

    public EnsembleController() throws IOException {
        int cores = Runtime.getRuntime().availableProcessors();
        System.out.println("Running " + System.getProperty("os.name") + " with " + cores + "CPUs");

        if (CLOptions.NUM_THREADS < 0) {
            cores = (cores > 2) ? cores / 2 : cores;
        } else {
            cores = CLOptions.NUM_THREADS;
        }

        System.out.println("Executor pool has " + cores + " workers");
        pool = Executors.newFixedThreadPool(cores);

        // assuming unique EOptions instanse is a unique calculation point, so no two
        // executors with the same options will be executed
        Set<EOptions> options = EOptions.readConfig(CONFIG_FILE);

        if (options.isEmpty()) {
            System.out.println("No valid options found in " + CONFIG_FILE);
            System.exit(2);
        }


        // setting main ensembles Lists
        for (EOptions opt : options) {
            IEnsemble ens = EnsemblesFactory.createEnsemble(opt);
            ensembles.add(ens);

            // fill mapping T -> ensembles list
            if (ensemblesResultValues.get(opt.getT()) == null) {
                ensemblesResultValues.put(opt.getT(), new LinkedHashMap<IEnsemble,
                        Deque<double[]>>());
            }

            ensemblesResultValues.get(opt.getT()).put(ens, new ArrayDeque<double[]>(3));
        }

        System.out.println("\n");
    }

    @Override
    public synchronized void stop() {
        notify();

        if (!running) return;
        running = false;

        System.out.println("Somebody stopping controller...");

    }

    public synchronized void start() throws InterruptedException {
        if (!running) return;

        System.out.println("Starting " + ensembles.size() + " ensembles' threads... ");
        ThreadLocalRandom rnd = ThreadLocalRandom.current();
        running = true;

        for (IEnsemble current : ensembles) {
            sleep(rnd.nextInt(5, 50)); // some start time distribution to avoid clutter
            // thread pool filled
            pool.execute(current);
            workerStates.put(current, true);
        }

        System.out.println("Ensembles started and running:");

        int i = 1;

        while (running) {
            drawStatus();

            // heavy status reports (saving energies & plots) only every 5th time
            if (i % 5 == 0) saveResults();

            refreshRunningStatus();
            i++;

            wait(CLOptions.REFRESH_DELAY);
        }


        for (IEnsemble current : ensembles) {
            current.stop();
        }

        sleep(500); // xxx: to be sure all threads stopped and to draw actual results on next line

        drawStatus();
        saveResults();

        System.out.println("Controller (" + ensembles.size() + " points) finished.");
        pool.shutdown();
        System.out.println("Thread pool shutted down.");
    }

    private void refreshRunningStatus() {
        // main ensembles refresh loop
        for (IEnsemble current : ensembles) {
            refreshStates(current);
        }
        stopIfFinished();
    }

    private void drawStatus() {
        System.out.println("\n\n"); // clear screen
        System.out.println(CLOptions.getOneLineSummary());
        System.out.println("-------------------------------------------------------------------------------------");

        // main ensembles refresh loop
        for (IEnsemble current : ensembles) {
            refreshResult(current);
            drawResult(current);
        }
    }

    /**
     * Save energy (for e ot gamma plot) for each temperature used
     */
    private void saveResults() {
        try {
            for (Integer currentKey : ensemblesResultValues.keySet()) {
                BufferedWriter writer = Files.newBufferedWriter(
                        GibbsConfigurationManager.getPath(currentKey + "K_" + CLOptions.NUM_PARTICLES + "pa_d" + CLOptions.MAX_DELTA_X + "_results" + ".txt"),
                        Charset.forName("UTF-8"));

                Set<IEnsemble> iEnsembles = ensemblesResultValues.get(currentKey).keySet();
                writer.write("#" + currentKey + "  " + iEnsembles);
                writer.newLine();

                for (IEnsemble ens : iEnsembles) {
                    final Deque<double[]> results = getResults(ens);
                    writer.write(getScientificResultString(results.peekLast()));
                    writer.newLine();
                }
                writer.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
            System.out.println("ERROR: can't save energy lists!");
        }
    }


    // =------------------ PRIVATE section -----------------
    private void stopIfFinished() {
        for (Boolean state : workerStates.values()) {
            if (state) return;
        }
        running = false;
        System.out.println("All threads seems finished work.");
    }

    private void refreshStates(IEnsemble current) {
        workerStates.put(current, !current.isFinished()); // true == running thread
    }

    private void drawResult(IEnsemble current) {
        boolean currentRunning = workerStates.get(current);
        String results = getSimpleResultsString(current);

        // warning – ensemble must implement toString() with meaningful tag for the next line
        // to write something useful
        System.out.println(current + "     \t#" + current.getCurrStep() +
                                   " (" + (int) (100 * ((float) current.getCurrStep() + 1) / current.getNumSteps())
                                   + "%)" + "\t[" + results + "]\t" + (currentRunning ? "" : " finished"));
    }


    private void refreshResult(IEnsemble current) {
        Deque<double[]> resValue = getResults(current); // we got list of energies for given ensemble
        final double[] results = current.getCurrentResult();

        if (!Arrays.equals(resValue.peekLast(), results)) {
            resValue.addLast(results);
        }
        if (resValue.size() > 3) resValue.pollFirst();
    }

    private String getScientificResultString(double[] result) {
        StringBuilder out = new StringBuilder();
        for (double value : result) {
            out.append(FULL_FORMAT.format(value));
            out.append("\t");
        }
        return out.toString();
    }

    private String getSimpleResultsString(IEnsemble current) {
        Deque<double[]> enrgs = getResults(current);
        StringBuilder out = new StringBuilder();

        for (double[] o : enrgs) {
            if (o.length > 2) {
                out.append(SHORT_FORMAT.format(o[0]));
                out.append("|");
                out.append(SHORT_FORMAT.format(o[1]));
                out.append(".");
            } else if (o.length > 1) {
                out.append(SHORT_FORMAT.format(o[0]));
                out.append("|");
                out.append(SHORT_FORMAT.format(o[1]));
            } else {
                out.append(SHORT_FORMAT.format(o[0]));
            }
            out.append(", ");
        }
        return out.toString();
    }

    private final Deque<double[]> getResults(IEnsemble current) {
        return ensemblesResultValues.get(current.getT()).get(current);
    }

}
