function [pbit0, pbit1] = llr2pb(llr)
pbit1 = exp(llr)./(1+exp(llr)); pbit0 = 1-pbit1;
