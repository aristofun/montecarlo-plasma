package com.butlitsky.mk;

/**
 * Pseudo potential ensemble. Electron-ion pseudopotential approximation
 * is unique for every [level numbers + temperature] combination.
 * <p/>
 * User: aristofun
 * Date: 02.03.13
 * Time: 17:19
 */
public class EnsemblePseudoPotential extends Ensemble {
//    private final PseudoDelegate myPseudo;

    public static final double LACING_POINT = 200.0;

    public EnsemblePseudoPotential(EOptions options) {
        super(options);
        initialize();
    }

    @Override
    /**
     * e^2 = 23,07077154753849 * e-20
     * Bohr (cm) = 5,2917721092e-9
     * k (in SGS) = 1,3806488 e-16
     *
     * => converted potential coeff. = e^2 / 7,30607881244045 * (e+5) = 3,
     * 15775016117465 e+5 = e+6/ 3,1668116505709 = 315775,01611746440408
     *
     * CURRENT version is for LEVEL 10-30 for 100K only
     * */
    protected final double getPotential(double r, boolean attraction) {

        if (attraction)   // ion-electron
        {
            if (r < LACING_POINT)
                return (-24.836465138387798 + 1.6042494381847763 * Math.pow(r, 0.3));
            else {
                return (-1 * 315775.01611746440408 / (T * r));
            }
        } else {
            if (r < 1)  // in Bor's radiuses
                return getPotential(1, false);
            else
                // The hard-coded Coloumb energy, always the same.
                return (315775.01611746440408 / (T * r));
        }
    }

    @Override
    protected double getPotentialAsym(double r, boolean ee, boolean ii) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    protected double getEnergy(double r, boolean attraction) {
        return getPotential(r, attraction);
    }

    @Override
    protected double getEnergyAsym(double r, boolean ee, boolean ii) {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

//    private abstract class PseudoDelegate {
//        abstract double pseudo(double x);
//    }
}
