import numpy as np
from connor_stevens_ode import connor_stevens_ode

def connor_stevens(t, I_ext):
    '''
    [V,X] = CONNOR_STEVENS(T, I_EXT) simulates connor_stevens neuron model
        the state variables are x = [V, n, m, h, a, b]
        the function returns
            1. V: Voltage trace
            2. x: state variables across time
    '''

    # Assumes that time is given in units of seconds, and converted into milisecond
    t_ms = 1000*t
    dt = t_ms[1] - t_ms[0]
    x = np.zeros(shape = (6, len(t_ms))) # V, n, m, h, a, b

    # initial conditions
    x[:, 0] = [-67.97, 0.1559, 0.01008, 0.9659, 0.5404, 0.2885]

    for tt in range(1, len(t_ms)):
        #print(tt, t_ms[tt])
        
        dx = connor_stevens_ode(x[:, tt-1], I_ext[tt])
        #print (dx)

        new_x = x[:, tt-1] + dt*dx;
        #print(new_x)

        # the voltage should be bounded between [-80, 50] mV
        new_V = new_x[0]
        if new_V > 50:
            new_V = 50
        elif new_V < -80:
            new_V = -80
        
        # the state variables should be bounded between [0, 1]
        new_x = new_x[1:]
        new_x[new_x<0] = 0
        new_x[new_x>1] = 1

        # Put the values back into x
        x[0, tt] = new_V
        x[1:, tt] = new_x

    V = x[0, :]
    x = x[1:, :]

    return V, x

if __name__ == "__main__": # run self-test
    t = np.arange(0, 0.2, 1e-5) # 200ms in 1ms steps
    I_ext = np.zeros_like(t)
    V, x = connor_stevens(t, I_ext)
    print(V)
    print(x)