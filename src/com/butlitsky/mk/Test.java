package com.butlitsky.mk;

import org.apache.commons.math3.util.FastMath;

import java.util.Arrays;

/**
 * Tests for different math stuff used in the programm etc.
 * <p/>
 * Date: 12.07.13
 * Time: 15:14
 * <p/>
 * Mac OS 10.9, Java 1.7, i7 winners:
 * Math.pow
 * FastMath.exp
 * Math.sqrt (tiny advantage over FastMath)
 * FastMath.cbrt
 * <p/>
 * Windows 8.1, Java 1.7, i5 winners:
 * Same results.
 */
public class Test {
    static double[][] array1, array2, tmp;

    public static final int DELAY = 507;

    public static void main(String[] srgs) throws InterruptedException {
        double f = -12E118;
        double q = f;
        // heatup
        for (int i = 0; i < 1000000; i++) {
            q += Math.pow(i * 0.0007, i * 0.000003);
            q -= 2 * FastMath.pow(i * 0.0007, i * 0.000003);
            q += StrictMath.pow(i * 0.0007, i * 0.000003);
            q += Math.exp(i * 0.0000007) - FastMath.exp(i * 0.0000007) + StrictMath.exp(i * 0.0000007);
        }
        array1 = new double[][]{
                {0.0, 1, 2, 3},
                {4, 5, 6, 7},
                {-1, -2, -3, -4}
        };

        array2 = new double[][]{
                {0.0, 0, 0, 0},
                {1, 1, 1, 1},
                {-2, -2, -2, -2}
        };

        System.out.println("array1:");
        outputArray(array1);
        System.out.println("array2:");
        outputArray(array2);

        tmp = array2;
        array2 = array1;
        array1 = tmp;
        tmp = null;

        System.out.println("\tafter switch");
        System.out.printf("tmp: " + tmp);
        System.out.println("array1:");
        outputArray(array1);
        System.out.println("array2:");
        outputArray(array2);

        /*
        Deque<Double>[] energies;

        energies = new Deque[2];
        energies[0] = new ArrayDeque<>(3);
        energies[1] = new ArrayDeque<>(5);

        energies[1].add(5.234);
        energies[1].add(-.234);
        energies[0].add(-.234);

        System.out.println("energies array: " + Arrays.toString(energies));
        System.out.println("0 and 1 samples: " + energies[0].getLast() + ", " + energies[1].getFirst());

        System.out.println("39/20 " + 39 / 20);
        System.out.println("20/20 " + 20 / 20);
        System.out.println("19/20 " + 19 / 20);
        System.out.println("39%20 " + 39 % 20);
        System.out.println("20%20 " + 20 % 20);
        System.out.println("19%20 " + 19 % 20);
        testPow();

        testDoubleExp();

        testLog();
        testDoubleSqrt();

        testCubicRoot(); */

        System.out.printf("\n\n");
    }

    private static void outputArray(double[][] array) {
        System.out.println(
                Arrays.toString(array[0]) + "\n" +
                        Arrays.toString(array[1]) + "\n" +
                        Arrays.toString(array[2]) + "\n"
        );
    }

    private static void testCubicRoot() throws InterruptedException {
        double f;
        long ms;
        final double one_third = (double) 1. / 3.;
        double initial = 1.29384023984E1;
        double coeff = 0.000102713;
        int steps = 70000000;

        System.out.println("\n\nTesting cubic root FastMath.cbrt");

// ------------------ StrictMath ---------------
        System.out.println("\n  StrictMath.cbrt: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += StrictMath.cbrt(i * coeff);
            f -= StrictMath.cbrt(++i * coeff);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

 /*       System.out.println("\n  StrictMath.pow(1/3): ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += StrictMath.pow(i * coeff, one_third);
            f -= StrictMath.pow(++i * coeff, one_third);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY); */

// ------------------ FastMath ---------------
        System.out.println("\n FastMath.cbrt: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += FastMath.cbrt(i * coeff);
            f -= FastMath.cbrt(++i * coeff);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

/*        System.out.println("\n  FastMath.pow(1/3): ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += FastMath.pow(i * coeff, one_third);
            f -= FastMath.pow(++i * coeff, one_third);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY); */

// -------------------- Math -----------------
        System.out.println("\n  Math.cbrt: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += Math.cbrt(i * coeff);
            f -= Math.cbrt(++i * coeff);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

/*        System.out.println("\n  Math.pow(1/3): ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += Math.pow(i * coeff, one_third);
            f -= Math.pow(++i * coeff, one_third);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);
  */


    }

    private static void testDoubleSqrt() throws InterruptedException {
        double f;
        long ms;
        double initial = 99.00717203801982302;
        double base = 0.000131793;
        int steps = 700000000;

        System.out.println("\n\nTesting double sqrt() ");

// ------------------ FastMath ---------------
        System.out.println("\n  FastMath: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += FastMath.sqrt(i * base);
            f -= FastMath.sqrt(++i * base);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

// ------------------ StrictMath ---------------
        System.out.println("\n  StrictMath: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += StrictMath.sqrt(i * base);
            f -= StrictMath.sqrt(++i * base);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

// -------------------- Math -----------------
        System.out.println("\n  Math: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += Math.sqrt(i * base);
            f -= Math.sqrt(++i * base);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);
    }

    private static void testDoubleExp() throws InterruptedException {
        double f;
        long ms;
        double initial = -1.62980238409242E3;
        double exp = 0.0000000717;
        int steps = 300000000;

        System.out.println("\n\nTesting double exp() ");

// ------------------ FastMath ---------------
        System.out.println("  FastMath: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += FastMath.exp(i * exp);
            f -= FastMath.exp(++i * exp);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

// ------------------ StrictMath ---------------
        System.out.println("\n  StrictMath: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += StrictMath.exp(i * exp);
            f -= StrictMath.exp(++i * exp);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

// -------------------- Math -----------------
        System.out.println("\n Math: ");

        ms = System.currentTimeMillis();
        f = initial;


        for (int i = 0; i < steps; i++) {
            f += Math.exp(i * exp);
            f -= Math.exp(++i * exp);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

    }

    private static void testPow() throws InterruptedException {
        double initial = -1.0200240928340982E2;
        double base = 0.00000197;
        double expo = 0.00000017;
        int steps = 70000000;


        System.out.println("Testing double pow() ");
        // -------------------- Math -----------------
        System.out.println("\n Math: ");
        double f = initial;

        long ms = System.currentTimeMillis();
        for (int i = 0; i < steps; i++) {
            f += Math.pow(i * base, i * expo);
            f -= Math.pow(++i * base, ++i * expo);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

// ------------------ FastMath ---------------
        System.out.println("\n  FastMath: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += FastMath.pow(i * base, i * expo);
            f -= FastMath.pow(++i * base, ++i * expo);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

// ------------------ StrictMath ---------------
        System.out.println("\n  StrictMath: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 0; i < steps; i++) {
            f += StrictMath.pow(i * base, i * expo);
            f -= StrictMath.pow(++i * base, ++i * expo);
        }

        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);
        Thread.sleep(DELAY);
    }

    private static void testLog() throws InterruptedException {
        double initial = -1.02002402;
        double base = 1.519702130123;
        int steps = 120000000;

        long ms;
        double f;

        System.out.println("Testing log() ");

        // -------------------- Math -----------------
        System.out.println("\n Math: ");
        f = initial;
        ms = System.nanoTime();

        for (int i = 1; i < steps; i++) {
            f += Math.log(i * base);
            f -= Math.log(++i * base);
        }
        System.out.print("" + (double) (System.nanoTime() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);

        // ------------------ FastMath ---------------
        System.out.println("\n  FastMath: ");

        ms = System.currentTimeMillis();
        f = initial;


        for (int i = 1; i < steps; i++) {
            f += FastMath.log(i * base);
            f -= FastMath.log(++i * base);
        }
        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);

        Thread.sleep(DELAY);


        // ------------------ StrictMath ---------------
        System.out.println("\n  StrictMath: ");

        ms = System.currentTimeMillis();
        f = initial;

        for (int i = 1; i < steps; i++) {
            f += StrictMath.log(i * base);
            f -= StrictMath.log(++i * base);
        }

        System.out.print("" + (double) (System.currentTimeMillis() - ms) / DELAY + "sec.\tf = " + f);
        Thread.sleep(DELAY);
    }
}
