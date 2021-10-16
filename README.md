# cstvCox
Codes for the paper:
Yanqing Sun,  Fei Heng (f.heng@unf.edu), Unkyung Lee, and Peter B. Gilbert, (2021) Estimation of conditional cumulative incidence functions under generalized semiparametric regression models with missing covariates, with application to biomarker correlates analyses in vaccine trials, Canadian Journal of Statistics

This folder contains the MATLAB code and a simulated data set:
* sim_data.mat: a simulated dataset
* MATLAB functions:
	+ main.m: main function
	+ esta.m, vara.m: compute estimates and estimated standard errors using two-stage AIPW method
	+ esti.m, vari.m: compute estimates and estimated standard errors using IPW method

In main.m, we use csvread() to load the simulated dataset: csvread('simdata.csv').

It includes input variables as shown below. 
	
INPUT:
------
	time: X=min{T,C}, where T is the failure time and C is the censoring time
	x: covariates with time-varying coefficients
	z1: covariates with time-constant coefficients
	z2: covariates with time-constant coefficients subject to missingness.
	rselec: selection indicator, 0 for z2 is missing, otherwise z2 is observed.
	A: auxialiry covariate
	cause: cause of failure V

	
PARAMETER SETTINGS:
-------------------
	dt: step of grid points
	times: grid points for t
	
OUTPUT:
-------
* Estimates of covariate coefficients: etahatall_a (AIPW), gamhat_a (AIPW), etahatall_i (IPW), gamhat_i (IPW)
* Estimated standard errors: sigmaeta_a (AIPW), sigmagam_a (AIPW), sigmaeta_i (IPW), sigmagam_i (IPW)
	
COMPUTATION TIME:
-----------------
The required computation time for sample size 800 is about 5 seconds running on the High Performance Computing Cluster at UNC Charlotte.
