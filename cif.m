function [F] = cif(t,lam,X,gam,Z) 
F = 1-exp(-((X*lam').*t).*exp(Z*gam'));
