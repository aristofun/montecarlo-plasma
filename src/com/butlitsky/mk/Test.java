package com.butlitsky.mk;

import org.apache.commons.math3.util.FastMath;

/**
 * Created with IntelliJ IDEA.
 * User: aristofun
 * Date: 12.07.13
 * Time: 15:14
 * To change this template use File | Settings | File Templates.
 */
public class Test {

    public static void main(String[] srgs) throws InterruptedException {
        System.out.println("Testing double pow() ");


        double f = - 12E118;
        double q = f;
        // heatup
        for (int i = 0; i < 1000000; i++) {
            q += Math.pow(i*0.0007, i*0.000003);
            q -= 2* FastMath.pow(i*0.0007, i*0.000003);
            q += StrictMath.pow(i*0.0007, i*0.000003);
            q += Math.exp(i*0.0000007) - FastMath.exp(i*0.0000007) + StrictMath.exp(i*0.0000007);
        }

// -------------------- Math -----------------
        System.out.println("  Math: ");

        long ms = System.currentTimeMillis();
        for (int i = 0; i < 10000000; i++) {
            f += Math.pow(i*0.0007, i*0.000003);
        }
        System.out.print("" + (double)(System.currentTimeMillis() - ms) / 1000 + "sec.  f = " + f);

        Thread.sleep(1000);

// ------------------ FastMath ---------------
        System.out.println("\n  FastMath: ");

        ms = System.currentTimeMillis();
        f = - 12E118;

        for (int i = 0; i < 10000000; i++) {
            f += FastMath.pow(i * 0.0007, i * 0.000003);
        }
        System.out.print("" + (double)(System.currentTimeMillis() - ms) / 1000 + "sec.  f = " + f);

        Thread.sleep(1000);

// ------------------ StrictMath ---------------
        System.out.println("\n  StrictMath: ");

        ms = System.currentTimeMillis();
        f = - 12E118;

        for (int i = 0; i < 10000000; i++) {
            f += StrictMath.pow(i * 0.0007, i * 0.000003);
        }
        System.out.print("" + (double)(System.currentTimeMillis() - ms) / 1000 + "sec.  f = " + f);

        Thread.sleep(1000);

        System.out.println("\nfinished\n");


        System.out.println("Testing double exp() ");

// -------------------- Math -----------------
        System.out.println("  Math: ");

        ms = System.currentTimeMillis();
        f = - 12E20;

        for (int i = 0; i < 100000000; i++) {
            f += Math.exp(i*0.0000007);
        }
        System.out.print("" + (double)(System.currentTimeMillis() - ms) / 1000 + "sec.  f = " + f);

        Thread.sleep(1000);

// ------------------ FastMath ---------------
        System.out.println("\n  FastMath: ");

        ms = System.currentTimeMillis();
        f = - 12E20;

        for (int i = 0; i < 100000000; i++) {
            f += FastMath.exp(i*0.0000007);
        }
        System.out.print("" + (double)(System.currentTimeMillis() - ms) / 1000 + "sec.  f = " + f);

        Thread.sleep(1000);

// ------------------ StrictMath ---------------
        System.out.println("\n  StrictMath: ");

        ms = System.currentTimeMillis();
        f = - 12E20;

        for (int i = 0; i < 100000000; i++) {
            f += StrictMath.exp(i*0.0000007);
        }
        System.out.print("" + (double)(System.currentTimeMillis() - ms) / 1000 + "sec.  f = " + f);

        Thread.sleep(1000);

        System.out.println("\nfinished");



        System.out.println("Testing double sqrt() ");

// -------------------- Math -----------------
        System.out.println("  Math: ");

        ms = System.currentTimeMillis();
        f = - 12E20;

        for (int i = 0; i < 100000000; i++) {
            f += Math.sqrt(i*0.00007);
        }
        System.out.print("" + (double)(System.currentTimeMillis() - ms) / 1000 + "sec.  f = " + f);

        Thread.sleep(1000);

// ------------------ FastMath ---------------
        System.out.println("\n  FastMath: ");

        ms = System.currentTimeMillis();
        f = - 12E20;

        for (int i = 0; i < 100000000; i++) {
            f += FastMath.sqrt(i*0.00007);
        }
        System.out.print("" + (double)(System.currentTimeMillis() - ms) / 1000 + "sec.  f = " + f);

        Thread.sleep(1000);

// ------------------ StrictMath ---------------
        System.out.println("\n  StrictMath: ");

        ms = System.currentTimeMillis();
        f = - 12E20;

        for (int i = 0; i < 100000000; i++) {
            f += StrictMath.sqrt(i*0.00007);
        }
        System.out.print("" + (double)(System.currentTimeMillis() - ms) / 1000 + "sec.  f = " + f);

        Thread.sleep(1000);

        System.out.println("\nfinished");
    }
}
