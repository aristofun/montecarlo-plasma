montecarlo-plasma
=================

Simple Markov chain Monte Carlo implementation for two component Coulomb plasma.
Periodical boundary conditions.
See more at: http://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo

This code uses good random algorithm: http://www.cs.gmu.edu/~sean/research/

It's made for multithreaded running of a set of points simultaneously.

Written for Java 7.

BTW, it works fast, really fast!

I mean, it's faster than the same code written in C++ and compiled by GCC 4, both on Mac OS X (Intel Core i7) and Windows 7 (Inter Core 2 Duo).

Contact me: aristofun at ya.ru


How to Run
===

		runmk.bat [OPTIONS]  – on Windows

		./runmk.comman [OPTIONS] – on Mac

Available options:

		-h							show this help and exit
 		-d,--delta <DELTA_FACTOR>	maxDx coeff. (1.3 default) for random particle shifts
 		-pa,--particles <NUM>		number of particles, if set all mk_config.ini options are ignored
		-po,--polka <POLKA>			polochka parameter value (2.0 default) / not available for custom potentials
		-r,--refresh <SECONDS>		console status refresh interval (20 sec. default)
		-w,--workers <NUM>			number of parallel points to calculate (default is MAX(2, CPUs/2))

Main input data file is mk_config.ini, it contains points to be calculated. Each line defines a point. They will be calculated in the same order, results are saved in subfolders.
Parameters must be delimited by any space character. Note that some global options above overrides individual points parameters if set.

Line format:

		T, Density, maxDx/y/z coeff., numParticles, numSteps, strategy (1 – save long tail configurations), resume existing config
