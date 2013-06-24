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
		
		maxDx/y/z 
			– if zero __and global --delta not defined__ used the whole box size (each trial move can place particle at any place inside the box).
			– if float number >= 1  __and global --delta not defined__ trial shift delta position == this number * average distance between particles (depends on density) along every axis.
			– if float number < 1 __and global --delta not defined__ trial shift delta is this fraction of a box size along each axis (i. e. 0.5 means trial particle shifts to half of box size along each axis). 
	
	
How it works
==

После распоковки в чистую директорию в ней можно запускать из командной строки.

Основной входной источник данных – файл mk_config.ini, он содержит набор точек для расчета, точки будут считаться по очереди в несколько параллельных потоков пока не пройдут указанное число шагов.

Каждые неск. секунд в консоль будет выводится текущий статус – на каком шаге находится каждый поток и последние пять значений удельной энергии этой конфигурации.

Кол-во потоков по умолчанию равно половине числа ядер в системе или двум (смотря что больше).

Программа довольно сильно греет процессор и заставляет на полную работать вентилятор.
Частично снизить эту проблему (например, чтобы на ночь поставить работать) можно задав в опциях запуска
параметр '-w 1', чтобы работал только один поток.

Ставить потоков больше числа логических ядер процессора нет смысла, общая скорость от этого не вырастет.

Набор конфигураций сохраняется в файле all_configs.dat каждой точки (в соотв. папке) через каждые 500 шагов.
Текущая конфигурация хранится в файле config.dat (первая строка там – Номер шага, уд. энергия, гамма,
далее – набор координат частиц).

Можно выключить прогу (по Ctrl + C, чтобы она gracefully закрылась), при след. запуске все точки, где указан последний параметр true будут продолжены с того места, где остановились (если не нарушится структура директорий где выполняется прога).
