package com.butlitsky.mk.ensembles;

import com.butlitsky.mk.options.EOptions;

import static java.lang.Math.pow;

/**
 * Pseudo potential ensemble. Electron-ion pseudopotential approximation
 * is unique for every [level numbers + temperature] combination.
 * <p/>
 * User: aristofun
 * Date: 02.03.13
 * Time: 17:19
 */
public class NVTEnsemblePseudoPotential extends NVTEnsemble {

    public static final double LACING_POINT = 185.0;

    public NVTEnsemblePseudoPotential(EOptions options) {
        super(options);

        System.out.println(" Pseudo Potential for T = " + options.getT() + ", " +
                "lacing pt. = " + LACING_POINT + ", created!");
    }

    /**
     * e^2 = 23,07077154753849 * e-20
     * Bohr (cm) = 5,2917721092e-9
     * k (in SGS) = 1,3806488 e-16
     * <p/>
     * U = e^2 / (r * Bohr * k)
     * <p/>
     * => converted potential coeff. = e^2 / 7,30607881244045 * (e+5) =
     * 3,15775016117465 e+5 = e+6/ 3,1668116505709 = 315775,01611746440408
     */
    protected final double getPotential(double r, boolean attraction) {

        if (attraction)   // ion-electron
        {
            if (r < LACING_POINT) {
                // 100 K, lvl 10-30 pseudo potential
                return -24.836465138387798 + 1.6042494381847763 * pow(r, 0.3);
            } else {
                return -1 * SCALE_FACTOR / (T * r);
            }
        } else {
            if (r < 2)  // in Bor's radiuses
                return getPotential(2, false);
            else
                // The hard-coded Coloumb energy, always the same.
                return (SCALE_FACTOR / (T * r));
        }
    }


    protected final double getEnergy(double r, boolean attraction) {
        // Assume non-coulomb potential makes zero contribution to Energy
        if (attraction && (r < LACING_POINT)) {
            return 0;
        } else {
            return getPotential(r, attraction); // xxx: temporary check full potential calculation
        }
    }


    @Override
    /**
     * used to introduce not symmetric i-i & e-e potentials
     *
     */
    protected double getPotentialAsym(double r, boolean ee, boolean ii) {
//  Barker approx.
//        if (ee)
//            return getPotential(r, false) * (1 - exp(-8.35E-4 * r * pow(T, 0.625)));
//
//        if (ii)
        if (ee || ii)
            return getPotential(r, false);
        else
            return getPotential(r, true);
    }


    protected final double getEnergyAsym(double r, boolean ee, boolean ii) {
// Barker approx.
//        if (ee)
//            return getPotential(r, false)
//                    * (1 - exp(-8.35E-4 * r * pow(T, 0.625)) * (1 - r * (8.35E-4) * pow(T, 0.625) * 0.625));
//
//        if (ii)
        if (ee || ii)
            return getEnergy(r, false);
        else
            return getEnergy(r, true);
    }
}
