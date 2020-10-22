function [s,cond] = NMDASynapse_step(dt,NT,Vpost,Mg,gmax,s)
E_1_2 = 16.13*log(Mg/3.57);
G_NMDA = (1+exp(-(Vpost-E_1_2)/16.13))^(-1);

alpha = 0.072; % mM-1ms-1
beta = 0.0066; % ms-1

dsdt = alpha*NT*(1-s)-beta*s;
s = s + dt*dsdt;

cond = G_NMDA*s*gmax;
    
return