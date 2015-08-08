package com.butlitsky.mk.ensembles;

import com.butlitsky.mk.options.CLOptions;
import com.butlitsky.mk.options.EOptions;
import org.apache.commons.math3.util.FastMath;

/**
 * Gibbse double box ensemble with Lennard-Johns potential
 * <p/>
 * myEpsilon HERE ACTS as Lennard-Johnes Epsilon energy factor!
 */
public class GibbsEnsembleLJ2 extends GibbsEnsemble {
    private static final double MINIMAL_DISTANCE = Double.MIN_VALUE * 100;
    private final double myEpsilon;

    /**
     * fixed initial sigma depending on density and rostar LJ parameter
     * in Bohrs radiuses!
     */
    private final double mySigma;

    protected GibbsEnsembleLJ2(EOptions options) {
        super(options,
              CLOptions.N_RESOLUTION_STEPS,
              CLOptions.N_RESOLUTION_STEPS * 3,
              CLOptions.N_RESOLUTION_STEPS * 5);

        myEpsilon = CLOptions.POLOCHKA;

//        ini file density is taken as ro1 reference parameter
        mySigma = FastMath.cbrt(CLOptions.RO_STAR1 / (options.getDensity())) / BOHR;

        System.out.println("\n T* = " + SHORT_FORMAT.format(1. / myEpsilon) +
                                   ", ro*1 = " + SHORT_FORMAT.format(CLOptions.RO_STAR1) +
                                   ", ro*2 = " + SHORT_FORMAT.format(CLOptions.RO_STAR2) +
                                   ", sigma = " + SHORT_FORMAT.format(mySigma));

        N = CLOptions.N1 + CLOptions.N2;
        Nei = N / 2;
        T = opt.getT();
        boxBorder = CLOptions.N1 / 2;

        System.out.println("\n N = " + N + ", Nei = " + Nei + ", N2 = " + CLOptions.N2 + ", boxBorder = " + boxBorder);

        // V is in CM^3 !!!
        V[0] = FastMath.pow(mySigma * BOHR, 3) * CLOptions.N1 / CLOptions.RO_STAR1;
        V[1] = FastMath.pow(mySigma * BOHR, 3) * CLOptions.N2 / CLOptions.RO_STAR2;
        Volume = V[0] + V[1];

        System.out.println("\n V[0] = " + V[0] + ", V[1] = " + V[1] + ", Volume = " + Volume);

        maxDXfactor = opt.getMaxDelta();
        maxDeltaV = CLOptions.MAX_DELTA_V;

        updateBoxSizes(false);
        updateLengths();
        updateDeltaX();

        System.out.println(myFolder + ": boxSize[0] = " + SHORT_FORMAT.format(boxSize[0]) +
                                   ", boxSize[1] = " + SHORT_FORMAT.format(boxSize[1]));

        System.out.println(
                " deltaX1 = " + SHORT_FORMAT.format(deltaX[0]) +
                        " deltaX2 = " + SHORT_FORMAT.format(deltaX[1])
        );

        initParticlesConfig();
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

    @Override
    protected void setStepType(int step) {
        // according to Athanatheos — initial steps to equilibrate boxes separately
        if (step < CLOptions.INITIAL_STEPS) {
            lastStepType = 0;
        } else {
            // otherwise x% of steps must be interchange particles only
            if (CLOptions.SWITCH_RATE > myRandom(1.0))
                lastStepType = 2;
            else
                lastStepType = ((step % 50) > 0) ? 0 : 1;
//                lastStepType = step % 2;
        }
    }
}