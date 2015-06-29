package com.butlitsky.mk.ensembles;

import com.butlitsky.mk.options.CLOptions;
import com.butlitsky.mk.options.EOptions;
import org.apache.commons.math3.util.FastMath;

/**
 * Gibbse double box ensemble with Lennard-Johns potential
 * <p/>
 * myEpsilon HERE ACTS as Lennard-Johnes Epsilon energy factor!
 */
public class GibbsEnsembleLJ extends GibbsEnsemble {
    private static final double MINIMAL_DISTANCE = Double.MIN_VALUE * 100;
    private final double myEpsilon;

    /**
     * fixed initial sigma depending on density and rostar LJ parameter
     * in Bohrs radiuses!
     */
    private final double mySigma;

    protected GibbsEnsembleLJ(EOptions options) {
        super(options);
        myEpsilon = CLOptions.POLOCHKA;
        mySigma = FastMath.cbrt(CLOptions.RO_STAR / (2 * options.getDensity())) / BOHR;

        System.out.print("\n T* = " + SHORT_FORMAT.format(1. / myEpsilon) + ", ro* = " + SHORT_FORMAT.format(CLOptions.RO_STAR)
                                 + ", sigma = " + SHORT_FORMAT.format(mySigma) + "\n");
    }

    /**
     * Pure Lennard–Johnes (r — in Bohr radiuses
     * Result must be dimensionless in kT values!
     */
    protected final double getPotential(double r) {
        if (r < MINIMAL_DISTANCE) {
            return getPotential(MINIMAL_DISTANCE);
        } else {
            final double sr = mySigma / r;
            return myEpsilon * 4.0 * (FastMath.pow(sr, 12) - FastMath.pow(sr, 6));
        }
    }

    /**
     * @param density - in cm-3
     * @return Ro* for Lennard–Johnes ensemble
     */
    @Override
    protected double scndryParam(double density) {
        return density * mySigma * mySigma * mySigma * BOHR * BOHR * BOHR; // pow(X, 3) optimization
    }

    /**
     * @param box — box number to check
     * @return Ro* for Lennard–Johnes ensemble
     */
    @Override
    protected double getSizeParam(int box) {
        return scndryParam(getDensitiesAvg()[box]);
    }

    @Override
    /**
     * used to introduce not symmetric i-i & e-e potentials
     *
     */
    protected double getPotentialAsym(double r, boolean ee, boolean ii) {
        return getPotential(r);
    }


    protected final double getEnergyAsym(double r, boolean ee, boolean ii) {
        return getEnergy(r);
    }


    protected final double getEnergy(double r) {
        return getPotential(r);
    }
}