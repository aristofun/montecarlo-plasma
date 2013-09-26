package com.butlitsky.mk;

/**
 * User: aristofun
 * Date: 02.03.13
 * Time: 17:19
 */
public class EnsemblePolochka extends Ensemble {
    /**
     * polochka deepness in kT
     */
    public static double EPSILON = 4.0;

    public static final double El = 4.80320427E-10;

    /**
     * scaling factor to convert energy value to kT units
     * SCALE_FACTOR == e^2 / (Bohr * k)
     */
    public static final double SCALE_FACTOR = 315775.01611746440408;

    protected final double myEpsilon;

    public EnsemblePolochka(EOptions options, boolean runinit) {
        super(options);
        myEpsilon = EPSILON;

        System.out.print(" polochka size = " + SHORT_FORMAT.format(SCALE_FACTOR / (T * myEpsilon)) + "\n");
        if (runinit) initialize();
    }

    @Override
    /**
     * e^2 = 23,07077154753849 * e-20
     * Bohr (cm) = 5,2917721092e-9
     * k (in SGS) = 1,3806488 e-16
     *
     * U = e^2 / (r * Bohr * k)
     *
     * => converted potential coeff. = e^2 / 7,30607881244045 * (e+5) =
     * 3,15775016117465 e+5 = e+6/ 3,1668116505709 = 315775,01611746440408
     * */
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
    protected final double getEnergy(double r, boolean attraction) {
        // – constant potential makes zero contribution to Energy
        if (attraction && (r < (SCALE_FACTOR / (T * myEpsilon)))) {
            return 0;
        } else {
            return getPotential(r, attraction);
        }
    }


}
