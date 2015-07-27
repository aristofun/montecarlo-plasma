package com.butlitsky.mk.ensembles;

import com.butlitsky.mk.options.CLOptions;
import com.butlitsky.mk.options.EOptions;

/**
 * Gibbse double box ensemble with polochka
 * <p/>
 * Created by aristofun on 25.09.14.
 */
public class GibbsEnsemblePolochka extends GibbsEnsemble {
    protected final double myEpsilon;

    protected GibbsEnsemblePolochka(EOptions options) {
        super(options);
        myEpsilon = CLOptions.POLOCHKA;

        System.out.print(" polochka=" + SHORT_FORMAT.format(SCALE_FACTOR / (T * myEpsilon)));
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
            if (r < (SCALE_FACTOR / (T * myEpsilon)))
                return (-1 * myEpsilon); //
            else {
                return (-1 * SCALE_FACTOR / (T * r));
            }
        } else {
            if (r < 1)  // in Bor's radiuses
                return getPotential(1, false);
            else
                // The hard-coded Coloumb energy, always the same.
                return (SCALE_FACTOR / (T * r));
        }
    }

    @Override
    protected double deltaTrialPotential() {
        return deltaTrialPotentialCoulomb();
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


    protected final double getEnergy(double r, boolean attraction) {
        // – constant potential makes zero contribution to Energy
        if (attraction && (r < (SCALE_FACTOR / (T * myEpsilon)))) {
            return 0;
        } else {
            return getPotential(r, attraction);
        }
    }
}