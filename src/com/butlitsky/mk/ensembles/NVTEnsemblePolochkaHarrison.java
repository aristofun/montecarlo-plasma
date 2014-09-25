package com.butlitsky.mk.ensembles;

import com.butlitsky.mk.options.CLOptions;
import com.butlitsky.mk.options.EOptions;

import static org.apache.commons.math3.util.FastMath.pow;

/**
 * Harrison R spheres summation summation algorithm.
 * <p/>
 * See page PHYSICAL REVIEW B 73, 212103, 2006
 * Simple calculation of Madelung constants, Walter A. Harrison
 * <p/>
 */
public class NVTEnsemblePolochkaHarrison extends NVTEnsemblePolochka {
    private final double myRcut;
    private final double myL;
    private final double myNcutoff2;
    private final int N;


    protected NVTEnsemblePolochkaHarrison(EOptions options) {
        super(options);

        N = CLOptions.HARRISON_N;
        myL = getBoxSize();
        myRcut = myL * (0.5 + N);
        myNcutoff2 = pow(Math.sqrt(2) + N, 2);
    }


    protected double getCurrentEnergy() {
        final int Nmax = N;
        final int N2 = N * N;
        final int NUM = getNumPart();
        final int NUM_2 = NUM / 2;
        final double Rcut2 = myRcut * myRcut;
        double result = 0;


        for (int i = 0; i < NUM; i++) {

//            double Q = 0;  // find effective Q for each particle

            // TODO: try effective i-e charges next (with double type)
            double energy = 0; // total energy of the i particle

            final int iQ = (i < NUM_2) ? 1 : -1; // FIRST NUM/2 are ions

            for (int nx = -Nmax; nx <= Nmax; nx++) {
                for (int ny = -Nmax; ny <= Nmax; ny++) {
                    for (int nz = -Nmax; nz <= Nmax; nz++) {
                        final double n2 = nx * nx + ny * ny + nz * nz;

                        // coarse grain sum limitation < R cutoff
                        if (n2 < myNcutoff2) {


                            // summ all other particles energy on i-th
                            for (int j = 0; j < NUM; j++) {
                                final double jx = Xs[j] + nx * myL;
                                final double jy = Ys[j] + nx * myL;
                                final double jz = Zs[j] + nx * myL;

                                // fine grain R cuttof (we count only charges inside R sphere)
                                if (n2 < N2 || (jx * jx + jy * jy + jz * jz < Rcut2)) {

                                    final int jQ = (j < NUM_2) ? 1 : -1;
                                    // exclude self i-i term in Zero Box
                                    if (i != j || nx != 0 || ny != 0 || nz != 0) {
                                        // i-j distance
                                        final double range = getRange(i, j, nx, ny, nz);
//                                        Q += getQ(iQ, jQ, range);
                                        energy += getEnergy(range, (iQ != jQ));
                                    }
                                }
                            }
                        }
                    }
                }

            }

//            final double corr = getCorrectionEnergy(iQ, Q);
//            if (corr != 0) {
//                System.out.println("energy: " + energy + ", Q: " + Q + ", iQ: " + iQ +
//                        "correction: " + corr);
//                System.out.println("Rcut " + myRcut);
//            }

//            energy += corr; // add compensating energy energy

            result += energy;
        }

//        System.out.println("RES: " + result);
        return result / 2.0;
    }

    /**
     * Compensating energy on the Rcutoff distance with given effective Koulobm charge
     *
     * @param iQ - particle for which calculating
     * @param q  – total effective charge (without probe particle!)
     * @return
     */
    private final double getCorrectionEnergy(int iQ, double q) {
        q += iQ;

        // TODO: check the sign on test simulations and returing value
        return 0.0 - iQ * q * getEnergy(myRcut, false);
    }

    /**
     * Get +1 or -1 or 0 effective charge of jq particle depending on the distance polochka
     * correction
     *
     * @param range
     * @return
     */
    private final double getQ(int iq, int jq, double range) {
        final int attract = iq * jq;
        final double polka = (SCALE_FACTOR / (T * myEpsilon));

        if ((attract < 0) && (range < polka)) {
            return jq * range / polka;
        } else {
            return jq;
        }

    }

    private final double getRange(int i, int j, int nx, int ny, int nz) {
        final double x = (Xs[j] - Xs[i]) + nx * myL;
        final double y = (Ys[j] - Ys[i]) + ny * myL;
        final double z = (Zs[j] - Zs[i]) + nz * myL;

        return Math.sqrt(x * x + y * y + z * z);
    }
}
