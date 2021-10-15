function [F] = cifest(lam,X,gam,Z) 
F = 1-exp(-(X*lam').*exp(Z*gam'));
