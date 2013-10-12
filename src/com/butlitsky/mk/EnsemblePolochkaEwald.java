package com.butlitsky.mk;

import static org.apache.commons.math3.special.Erf.erfc;
import static org.apache.commons.math3.util.FastMath.*;

/**
 * General Ewald summation algorithm.
 * <p/>
 * See page 74 in
 * Algorithms for Molecular Dynamics Simulations (Advancing the Computational Horizon)
 * by FREDRIK HEDMAN
 * <p/>
 * <a href="http://www.google.ru/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&ved=0CC0QFjAA&url=http%3A%2F%2Fwww.mmk.su.se%2Fdocuments%2Fpublications%2Fthesis_fredrik.pdf&ei=4UQsUsKRKeik4ATX_4CYBg&usg=AFQjCNHfS25g4lwiSOrjjcqHcQJMvST7GQ&sig2=dQn4yb3NGhebXCn_QAjWhQ&bvm=bv.51773540,d.bGE"
 * >link</a>
 */
public class EnsemblePolochkaEwald extends EnsemblePolochka {
    /**
     * Ewald alpha
     */
    public static double DELTA = 0.00000001;

    /**
     * Ewal n cutoff
     */
    public static int Ncutoff = 3;


    private static boolean notPrinted = true;

    private final double myDelta;
    private final int myNcutoff;
    private final double myAlpha;
    private final double myRcut;
    private final double myUself;


    public EnsemblePolochkaEwald(EOptions options) {
        super(options,false);
        myDelta = DELTA;
        myNcutoff = Ncutoff;
        myAlpha = PI * myNcutoff / (getBoxSize() * sqrt(-log(myDelta)));

        myRcut = PI * Ncutoff / (myAlpha * myAlpha * getBoxSize()); // getBoxSize();
        myUself = Uself();

        printStat(myDelta, myNcutoff, myAlpha, myAlpha * getBoxSize(), myRcut / getBoxSize());
        initialize();
    }

    private static void printStat(double delta, int ncutoff, double alpha, double alphaL, double rcut) {
        if (notPrinted) {
            notPrinted = false;
            System.out.print("EWALD delta = " + delta + ", Ncut = " + ncutoff);
            System.out.print(" EWALD alphaL <= " + alphaL + " EWALD alpha <= " + alpha);
            System.out.println("Rcut = " + rcut);
        }
    }

    protected double getCurrentEnergy() {
//        return super.getCurrentEnergy();
//        System.out.println("writing fixed energy!");
        final double ureal = Ureal();
        final double urecip = Urecip();
        final double ubc = Ubc();

//        System.out.println("Real: " + SHORT_FORMAT.format(ureal)
//                + ", Recip: " + SHORT_FORMAT.format(urecip)
//                + ", Myself: " + SHORT_FORMAT.format(myUself)
//                + ", BC: " + SHORT_FORMAT.format(ubc));

        return ureal
                + urecip
                - myUself
                + ubc
                ; // XXX try to turn off some parts
    }

    private final double Ureal() {
        double result = 0.;
        final double L = getBoxSize(); // all distances are in BOHR radiuses!!!
        final double alpha = myAlpha; // this Alpha is in 1/BOHR radiuses!
        final int Nmax = myNcutoff;
        final double rcut = myRcut;
        final int NUM = getNumPart();

        for (int Nx = -Nmax; Nx <= Nmax; Nx++) {
            for (int Ny = -Nmax; Ny <= Nmax; Ny++) {
                for (int Nz = -Nmax; Nz <= Nmax; Nz++) {

                    // go over particles
                    for (int i = 0; i < NUM; i++) {
                        for (int j = 0; j < NUM; j++) {
                            // exclude self term
                            if (Nx == 0 && Ny == 0 && Nz == 0 && i == j) continue;

                            final double range = getRange(i, j, Nx, Ny, Nz, L); // in BOHR radiuses

                            if (range < rcut) {
                                final double potential;

                                if (j < (NUM / 2))   // First _num/2 are IONS
                                {
                                    if (i < (NUM / 2)) // ION-ION
                                        potential = getEnergy(range, false);
                                    else              // ION - electron
                                        potential = getEnergy(range, true);
                                } else                   // Last _num/2 are Electrons
                                {
                                    if (i < (NUM / 2)) // Electron - ION
                                        potential = getEnergy(range, true);
                                    else              // Electron - Electron
                                        potential = getEnergy(range, false);
                                }
                                // (2.38) formula in book
                                result += potential * erfc(alpha * range);
                            }
                        }
                    }
                }
            }
        }

        return result / 2.0;
    }

    private final double getRange(int i, int j, int nx, int ny, int nz, double L) {
        final double x = fit2box(Xs[j] - Xs[i]) + nx * L;
        final double y = fit2box(Ys[j] - Ys[i]) + ny * L;
        final double z = fit2box(Zs[j] - Zs[i]) + nz * L;

        return StrictMath.sqrt(x * x + y * y + z * z);
    }

    private final double Urecip() {
        double result = 0.;
        final double L = getBoxSize(); // all distances are in BOHR radiuses!!!
        final double alphaL2 = myAlpha * L * myAlpha * L; // this AlphaL is UNITLESS
        final double Pi2 = PI * PI;
        final double Pi2_L = PI * 2 / L;
        final int Nmax = myNcutoff;
        final int Nmax2 = myNcutoff * myNcutoff;
        final double q = e;
        final int NUM = getNumPart();

        int n2 = 0;

        for (int Nx = -Nmax; Nx <= Nmax; Nx++) {
            for (int Ny = -Nmax; Ny <= Nmax; Ny++) {
                for (int Nz = -Nmax; Nz <= Nmax; Nz++) {
                    n2 = Nx * Nx + Ny * Ny + Nz * Nz;
                    if (n2 == 0) continue;

                    if (n2 < Nmax2) { // XXX check if this influence results
                        final double coeff = exp(-Pi2 * n2 / alphaL2) / n2;
                        double Qcos = 0.0, Qsin = 0.0;

                        for (int i = 0; i < NUM; i++) {
                            if (i < (NUM / 2)) {   // First _num/2 are IONS
                                Qcos += q * cos(Pi2_L * (Nx * Xs[i] + Ny * Ys[i] + Nz * Zs[i]));
                                Qsin += q * sin(Pi2_L * (Nx * Xs[i] + Ny * Ys[i] + Nz * Zs[i]));
                            } else {                  // Last _num/2 are Electrons
                                Qcos -= q * cos(Pi2_L * (Nx * Xs[i] + Ny * Ys[i] + Nz * Zs[i]));
                                Qsin -= q * sin(Pi2_L * (Nx * Xs[i] + Ny * Ys[i] + Nz * Zs[i]));
                            }
                        }
                        // (2.40) formula in book
                        result += coeff * (Qcos * Qcos + Qsin * Qsin);
                    }
                }
            }
        }

        return result / (2.0 * PI * L * k * T * BOHR);//XXX BOHR radiuses and kT?
    }

    private final double Uself() {
        return SCALE_FACTOR * getNumPart() * myAlpha / (sqrt(PI) * T);
    }

    private final double Ubc() {
        double resultx = 0.;
        double resulty = 0.;
        double resultz = 0.;
        final double Pi2_3L3 = PI * 2 / (3 * getBoxSize() * getBoxSize() * getBoxSize());
        // all distances are in BOHR radiuses!!!
        final double q = e;
        final int NUM = getNumPart();

        for (int i = 0; i < NUM; i++) {
            if (i < (NUM / 2)) {   // First _num/2 are IONS
                resultx += Xs[i];
                resulty += Ys[i];
                resultz += Zs[i];
            } else {
                resultx -= Xs[i];
                resulty -= Ys[i];
                resultz -= Zs[i];
            }
        }

        return SCALE_FACTOR
                * Pi2_3L3 * (resultx * resultx + resulty * resulty + resultz * resultz) / T;
    }

/*    private final double getEwPotential(int i, int j, double distance) {
        final double r = StrictMath.sqrt(dSquared(Xs[j] - Xs[i], Ys[j] - Ys[i], Zs[j] - Zs[i]));

        if (i != j) {
            if (j < (numPart / 2))   // First _num/2 are IONS
            {
                if (i < (numPart / 2)) // ION-ION
                    return getEnergyAsym(r, false);
                else              // ION - electron
                    return getEnergyAsym(r, true);
            } else                   // Last _num/2 are Electrons
            {
                if (i < (numPart / 2)) // Electron - ION
                    return getEnergyAsym(r, true);
                else              // Electron - Electron
                    return getEnergyAsym(r, false);
            }
        }
        return 0;
    }
    */
}
