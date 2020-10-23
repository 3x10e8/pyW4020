import numpy as np
from numpy import exp

# Define Default Hodgkin-Huxley Neuron ODEs. This function handle is 
    # written by Yiyin Zhou.
def hhn_ode(x, I_osc = 0):
    v = x[0]

    # compute dn
    a = exp(-(v+55)/10) -1
    dn = 1 -x[1]
    if a == 0:
        dn *= 0.1
    else:
        dn *= -0.01*(v+55)/a 
    dn -= x[1]*(0.125*exp(-(v+65)/80))

    # compute dm
    a = exp(-(v+40)/10) -1
    dm = 1 -x[2]
    if a == 0:
        dm *= 1 # no change
    else:
        dm *= -0.1*(v+40)/a
    dm -= x[2] * (4*exp(-(v+65)/18))

    dh = (1-x[3]) * (0.07*exp(-(v+65)/20))
    dh -= x[3] / (exp(-(v+35)/10)+1)

    dV = I_osc 
    dV -= 120 * (x[2]**3) * x[3] * (v-50) 
    dV -=  36 * (x[1]**4) * (v+77) 
    dV -= 0.3 * (v+54.387)

    dx = np.array([dV, dn, dm, dh])

    return dx

# TO DO, figure out matrix multiplication
# Define Default Rinzel Neuron model ODEs.
'''
def rinzel_ode(x):
    v = x[0]

    a = exp(-(v+55)/10)-1
    alpha_n = np.matmul([a==0 a!=0], [0.01; -0.01*(v+55)/a])

    a = exp(-(v+40)/10)-1;
    alpha_m = ([a==0 a~=0])*[1; -0.1*(v+40)/a]

    alpha_h = 0.07*exp(-(v+65)/20) ;                     

    beta_n = 0.125*exp(-(v+65)/80); 
    beta_m = 4*exp(-(v+65)/18);
    beta_h = 1/(exp(-(v+35)/10)+1);                  

    n_infty = alpha_n/(alpha_n + beta_n);
    m_infty = alpha_m/(alpha_m + beta_m);    
    h_infty = alpha_h/(alpha_h + beta_h);

    # compute dR
    S = 1.2714;  % s = (1-h_inf(0)) / n_inf(0) 
    R_infty = S/(1+S^2)*(n_infty + S*(1-h_infty));
    tau_R = 1 + 5*exp(-(v+55)^2/55^2);   
    dR = -3*(x(2)-R_infty) / tau_R;

    # Update the ionic currents and membrane voltage:
    dV = I_osc -120*m_infty^3*(1-x(2))*(v-50) - 36*(x(2)/S).^4*(v+77) - 0.3*(v+54.387);
    dx = [dV dR];

    return dx
'''

# Define Default Wilson Neuron model ODEs.
def wilson_ode(x):
    v = x[0]
    R_infty = 0.0135*v + 1.03
    dR = -(x[1] -R_infty) / 1.9
    dV = I_osc 
    dV -= (17.81 + 0.4771*v + 3.263e-3*v**2)*(v-55) 
    dV -= 26*x[1]*(v+92)

    dx = np.array([dV, dR])
    return dx

def hodgkin_huxley_ode(t, I_ext, Model='HH', State=0,
                        Coupling='Additive', ConstI=0, FullOutput=False):
    '''
    HODGKIN_HUXLEY_ODE simulate a given neuron model
       [out] = HODGKIN_HUXLEY_ODE(T,I_EXT,VARARGIN) simulates a given neuron 
       model with current input I_EXT for time course T. The neuron model is 
       described as a set of ordinary differential equations (ODEs). The
       resultant voltage trace of the neuron is return in OUT.
     
       [...] = HODGKIN_HUXLEY_ODE(...,'Model',MODEL) finds PRC for the given
       neuron model MODEL. MODEL can be one of 'HH', as default, 'Rinzel', and
       'Wilson', or a Matlab function handle for user defined model. When user 
       defined model is provided, please also use 
       HODGKIN_HUXLEY_PRC(...,'State',state) to specify the initial values of
       the state variables for the given model.
    
       [...] = HODGKIN_HUXLEY_ODE(...,'Coupling',COUPLING) simulates the 
       neuron under specified coupling method. When COUPLING is set to 
       'Additive' as default, the input current I_EXT is additively coupled to 
       the state variables of the neuron. On the other hand, when COUPLING is 
       set to 'MULTIPLICATIVE', the input current I_EXT is multiplicatively
       coupled to the state variables of the neuron.
       
       [...] = HODGKIN_HUXLEY_ODE(...,'CONSTI',CONST_I) injects a constant
       current bias CONST_I in addition to I_EXT into the neuron.
    
       [...] = HODGKIN_HUXLEY_ODE(...,'FULLOUTPUT', TRUE) returns not only the
       voltage trace but also the traces of all state variables in a N-by-M
       matrix out where N is the length of time course T, and M is the number
       of state variables.
    
       Author:      Chung-Heng Yeh <chyeh@ee.columbia.edu>
    
       Copyright 2010-2012   Chung-Heng Yeh
    '''

    # Handle the optional input parameters.
    # =====================================================================
    # Model: Specify the Neuron ODEs.
    # State: Specify the initial values of neuron state variables.
    # Coupling: Specify the coupling method.
    # FullOutput: Specify full output.
    # ConstI: Specify the constant bias.

    # Convert string inputs to lower case for easier handling
    Model = Model.lower()
    Coupling = Coupling.lower()

    # Handle the unexpected input parameter.
    if not Model in ['hh', 'wilson', 'rinzel']:
        print("Model can be either 'HH', 'Wilson', or 'Rinzel'.")
    if not isinstance(State, (int, float)):
        print('State must be an int or float.')
    if not Coupling in ['additive','multiplicative']:
        print("Coupling can be 'Additive' or 'Multiplicative'.")
    if not isinstance(FullOutput, bool):
        print('FullOutput must be a bool.')
    if not isinstance(ConstI, (int, float)):
        print('ConstI must be an int or float.')
    
    # Setup the numerical method.
    # =====================================================================
    # Initialize the Neuron function handel and state variables. 
    if Model.lower() == 'hh':       # Hodgkin-Huxley Model
        neu_model = hhn_ode
        neu_state = [0., 0., 0., 1.]
    elif Model.lower() == 'rinzel':  # Rinzel Model
        neu_model = rinzel_ode
        neu_state = [-60., 0.4]       
    elif Model.lower() == 'wilson': # Wilson Model
        neu_model = wilson_ode
        neu_state = [-70., -0.088]
    else:                           # User specified Model
        print ('Please see code for settig up a new model.')
        #neu_model = # provide an ode function
        #neu_state = # provide a state vector   

    # Make a numpy array for vector operations 
    neu_state = np.array(neu_state)

    # Assume that the time is given in seconds, convert time to millisecond.
    dt = 1000.0*(t[1]-t[0])

    # Initialize membrane voltage vector.
    out = np.zeros(shape=(len(t),len(neu_state)))

    # Simulate the neuron model.
    for n in range(len(I_ext)):
        # Compute the derivatives of neuron state variables.
        if neu_model == hhn_ode:
            # provide the extra contant I bias
            d_neu_state = neu_model(neu_state, I_osc = ConstI)
        else:
            d_neu_state = neu_model(neu_state)

        # Use different coupling method. 
        if Coupling == 'additive': # Additive Coupling
            d_neu_state[0] += I_ext[n]
        else: # Multiplicative Coupling
            d_neu_state *= I_ext[n]

        # Use the forward Euler method to integrate.
        neu_state += dt*d_neu_state

        # Export the membrane voltage to output.
        out[n,:] = neu_state

    if not FullOutput:
        out = out[:, 0]

    return out

if __name__ == '__main__':
    t = np.arange(0, 0.2, 1e-6)
    I_ext = 1e-6 * np.sin(t)
    Model = 'HH'
    out = hodgkin_huxley_ode(t, I_ext, Model='HH')
    print(out)