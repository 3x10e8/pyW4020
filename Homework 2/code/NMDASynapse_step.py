from numpy import log, exp

def NMDASynapse_step(dt,NT,Vpost,Mg,gmax,s):
	E_1_2 = 16.13*log(Mg/3.57) # this is natural log
	G_NMDA = (1+exp(-(Vpost-E_1_2)/16.13))**(-1)

	alpha = 0.072 # mM-1ms-1
	beta = 0.0066 # ms-1

	dsdt = alpha*NT*(1-s)-beta*s
	s += dt*dsdt
	cond = G_NMDA*s*gmax
	    
	return s, cond

if __name__ == '__main__':
	dt = 1e-6
	NT = 1
	Vpost = 0
	Mg = 1
	gmax = 1
	s = 0
	s, cond = NMDASynapse_step(dt,NT,Vpost,Mg,gmax,s)
	print(s, cond)