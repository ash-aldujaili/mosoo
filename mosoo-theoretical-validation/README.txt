##################################
# MO-SOO Theoretical Validation  #
##################################

This directory contains scripts/data that are used
in validating the theoretical bounds of MO-SOO with respect 
to eight instances of a bi-objective problem (n = 1, 2; Psi take one of four values). These instances have been encoded with n%d-i%d to denote the ith problem in the nth dimension.

## In essence:
	1. computeBoundMOSOO.m: computes the theoretical bounds for the loss and indicator using Symbolic Maths
	2. compileTheoreticalPlots.py: a python scripts that use all the data generated and other scripts(*.py, *.so files) to compile a set of tables that
	   generated the theoretical validation plots in the journal
	3. pf_data/pf-n%d-mosoo-i%d.txt (dimension, instance): contains numerically-obtained pareto fronts
	4. pf_data/ext-front-n%d-mosoo-i%d (dimension, instance): contains the extrema of the pareto front-n%d-mosoo-i%d
	5. numerical_data/resultfile-f1-d%d-r%d-fx.txt (dimension, instance): contains the numerical simulation results of MO-SOO on the eight instances where each row represents
		the tuple (timestamp, objective1-value, objective2-value)  only those samplings that dominated the previous points are recorded here
	6. psi_data/psi-n%d-mosoo.txt: Psi (conflict dimensionality) values for each instance/ dimension (can be computed by compileTheoreticalPlots.py)
	7. numerical_data/ind-n%d-mosoo-i%d : numerically-computed indicator values
	8. numerical_data/loss-n%d-mosoo-i%d : numerically-computed loss values
	9. theoretical_data/ 
	10. hmax_data/hmax-n%d-mosoo-i%d : numerically-computed hmax values for MO-SOO on different instances of the problem

|=/\=| P.S: *.so files are called in python with or without their extension depending on the platform. You may want to remove/insert .so in the lines specifying
		the library files in paretofront.py and epsilonindicator.py

## To generate the tables from the available data:
>> python compileTheroticalPlots.py

## To generate the numerical results again and generate the tables
>> matlab run_me.m

## To generate the theoretical results again
>> matlab computeBoundMOSOO.m
	

Tested on Windows 7, with MATLAB 2015b.

Abdullah Al-Dujaili, 2016