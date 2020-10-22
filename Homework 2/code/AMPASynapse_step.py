def AMPASynapse_step(dt,NT,gmax,s):
	alpha = 1.1 # mM-1ms-1
	beta = 0.19 # ms-1

	dsdt = alpha*NT*(1-s) -beta*s
	s += dt*dsdt
	cond = s*gmax
	    
	return s, cond


if __name__ == '__main__': # run self-test
	dt = 1e-6
	NT = 1
	gmax = 1e-3
	s = 0
	s, cond = AMPASynapse_step(dt,NT,gmax,s)
	print ('state: {}'.format(s))
	print ('conductance: {}'.format(cond))