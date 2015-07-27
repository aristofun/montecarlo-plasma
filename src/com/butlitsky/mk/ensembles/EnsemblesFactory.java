package com.butlitsky.mk.ensembles;

import com.butlitsky.mk.IEnsemble;
import com.butlitsky.mk.options.CLOptions;
import com.butlitsky.mk.options.EOptions;

/**
 * Ensembles only created on the Factory to ensure state is correctly loaded before running.
 * <p/>
 * Created by aristofun on 21.07.14.
 */
public class EnsemblesFactory {
    private EnsemblesFactory() {
    }

    public static IEnsemble createEnsemble(EOptions opt) {
        IEnsemble ensemble;
        switch (CLOptions.ENSEMBLE_TYPE) {
            case 1:
                ensemble = new NVTEnsemblePolochkaEwald(opt);
                ensemble.loadState();
                break;
            case 2:
                ensemble = new NVTEnsemblePolochkaHarrison(opt);
                ensemble.loadState();
                break;
            case 3:
                ensemble = new NVTEnsemblePseudoPotential(opt);
                ensemble.loadState();
                break;
            case 4:
                ensemble = new GibbsEnsemblePolochka(opt);
                ensemble.loadState();
                break;
            case 5:
                ensemble = new GibbsEnsembleLJ(opt);
                ensemble.loadState();
                break;
            case 6:
                ensemble = new GibbsEnsembleLJ2(opt);
                ensemble.loadState();
                break;
            default:
                ensemble = new NVTEnsemblePolochka(opt);
                ensemble.loadState();
                break;
        }
        return ensemble;
    }
}
