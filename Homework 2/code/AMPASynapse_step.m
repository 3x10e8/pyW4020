function [s,cond] = AMPASynapse_step(dt,NT,gmax,s)
alpha = 1.1; % mM-1ms-1
beta = 0.19; % ms-1

dsdt = alpha*NT*(1-s)-beta*s;
s = s + dt*dsdt;

cond = s*gmax;
    
return