package com.butlitsky.mk;

/**
 * Created with IntelliJ IDEA.
 * User: aristofun
 * Date: 02.03.13
 * Time: 17:08
 */
public interface IEnsembleController {

    /** starts current config to execute */
    void start() throws InterruptedException;

    /** gracefully stop all controller threads and return then */
    void stop();
}
