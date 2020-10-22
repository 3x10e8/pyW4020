import numpy as np
from numpy import exp

def connor_stevens_ode(x, I_ext):
    '''
    [DX] = CONNOR_STEVENS_ODE(X, I_EXT) 
    Computes the gradient of all states given 
        previous state value X and the input I_EXT
    The state variables are x = [V, n, m, h, a, b]
    '''

    ms=-5.3;  ns=-4.3; hs=-12;
    gNa=120; gK=20;  gL=0.3; ga=47.7;
    ENa=55;  EK=-72; EL=-17; Ea=-75;
    
    v = x[0];
    n = x[1]; m = x[2]; h = x[3];
    a = x[4]; b = x[5];
    
    alpha = exp(-(v+50+ns)/10)-1
    if abs(alpha) <= 1e-7:
        alpha = 0.1
    else:
        alpha = -0.01*(v+50+ns)/alpha

    beta = .125*exp(-(v+60+ns)/80)
    n_inf = alpha/(alpha+beta)
    tau_n = 2/(3.8*(alpha+beta))

    alpha = exp(-(v+35+ms)/10)-1.
    if abs(alpha) <= 1e-7:
        alpha = 1
    else:
        alpha = -.1*(v+35+ms)/alpha

    beta = 4*exp(-(v+60+ms)/18)
    m_inf = alpha/(alpha+beta)
    tau_m = 1/(3.8*(alpha+beta))

    alpha = .07*exp(-(v+60+hs)/20)
    beta = 1/(1+exp(-(v+30+hs)/10))
    h_inf = alpha/(alpha+beta)
    tau_h = 1/(3.8*(alpha+beta))

    a_inf = (.0761*exp((v+94.22)/31.84)/(1+exp((v+1.17)/28.93)))**(1/3)
    tau_a = .3632+1.158/(1+exp((v+55.96)/20.12))
    b_inf = (1/(1+exp((v+53.3)/14.54)))**4
    tau_b = 1.24+2.678/(1+exp((v+50)/16.027))

    i_na = gNa * (m**3) * h * (v - ENa)
    i_k = gK * (n**4) * (v - EK)
    i_l = gL * (v - EL)
    i_a = ga * (a**3) * b * (v - Ea)

    d_v = I_ext - i_na - i_k - i_l - i_a
    d_n = (n_inf-n)/tau_n
    d_m = (m_inf-m)/tau_m
    d_h = (h_inf-h)/tau_h
    d_a = (a_inf-a)/tau_a
    d_b = (b_inf-b)/tau_b

    dx = np.array([d_v, d_n, d_m, d_h, d_a, d_b])

    return dx

if __name__ == "__main__": # run self-test
    V=1; n=1; m=1; h=1; a=1; b=1
    x = [V, n, m, h, a, b]
    I_ext = 0
    dx = connor_stevens_ode(x, I_ext)
    print(dx)