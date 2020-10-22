%HODGKIN_HUXLEY_ODE simulate a given neuron model
%   [out] = HODGKIN_HUXLEY_ODE(T,I_EXT,VARARGIN) simulates a given neuron 
%   model with current input I_EXT for time course T. The neuron model is 
%   described as a set of ordinary differential equations (ODEs). The
%   resultant voltage trace of the neuron is return in OUT.
% 
%   [...] = HODGKIN_HUXLEY_ODE(...,'Model',MODEL) finds PRC for the given
%   neuron model MODEL. MODEL can be one of 'HH', as default, 'Rinzel', and
%   'Wilson', or a Matlab function handle for user defined model. When user 
%   defined model is provided, please also use 
%   HODGKIN_HUXLEY_PRC(...,'State',state) to specify the initial values of
%   the state variables for the given model.
%
%   [...] = HODGKIN_HUXLEY_ODE(...,'Coupling',COUPLING) simulates the 
%   neuron under specified coupling method. When COUPLING is set to 
%   'Additive' as default, the input current I_EXT is additively coupled to 
%   the state variables of the neuron. On the other hand, when COUPLING is 
%   set to 'MULTIPLICATIVE', the input current I_EXT is multiplicatively
%   coupled to the state variables of the neuron.
%   
%   [...] = HODGKIN_HUXLEY_ODE(...,'CONSTI',CONST_I) injects a constant
%   current bias CONST_I in addition to I_EXT into the neuron.
%
%   [...] = HODGKIN_HUXLEY_ODE(...,'FULLOUTPUT', TRUE) returns not only the
%   voltage trace but also the traces of all state variables in a N-by-M
%   matrix out where N is the length of time course T, and M is the number
%   of state variables.
%
%   Author:      Chung-Heng Yeh <chyeh@ee.columbia.edu>
%
%   Copyright 2010-2012   Chung-Heng Yeh

function out = hodgkin_huxley_ode(t,I_ext,varargin)    
    % Handle the optional input parameters.
    % =====================================================================
    p = inputParser;
    p.KeepUnmatched = true;
    % Specify the Neuron ODEs.
    addParamValue(p,'Model','HH',@(x) isa(x,'function_handle') || any(validatestring(x,{'HH','Wilson','Rinzel'})));
    % Specify the initial values of neuron state variables.
    addParamValue(p,'State',0, @isnumeric);
    % Specify the coupling method.
    addParamValue(p,'Coupling','Additive',@(x) any(validatestring(x,{'Additive','Multiplicative'})));
    % Specify full output.
    addParamValue(p,'FullOutput',false,@islogical);
    % Specify the constant bias.
    addParamValue(p,'ConstI',0,@isnumeric);
    
    % Parse the optional parameters.
    parse(p,varargin{:});

    % Handle the unexpected input parameter.
    UnmatchedParam = fieldnames(p.Unmatched);
    if ~isempty(UnmatchedParam)
        error(['"',UnmatchedParam{1},'" is not a valid parameter.']);
    end
    
    % Setup the numerical method.
    % =====================================================================
    % Initialize the Neuron function handel and state variables. 
    if strcmpi( p.Results.Model, 'HH' ),           % Hodgkin-Huxley Model
        neu_model = @hhn_ode;
        neu_state = [0 0 0 1];
    elseif strcmpi( p.Results.Model, 'Rinzel' ),   % Rinzel Model
        neu_model = @rinzel_ode;
        neu_state = [-60 0.4];        
    elseif strcmpi( p.Results.Model, 'Wilson' ),   % Wilson Model
        neu_model = @wilson_ode;
        neu_state = [-70 -0.088];        
    else                                           % User specified Model
        neu_model = p.Results.Model;
        neu_state = p.Results.State;         
    end
    % Use different coupling method. 
    if strcmpi( p.Results.Coupling , 'Additive'), % Additive Coupling
        neu_func = @(state, I) neu_model(state) + [I, zeros(1,numel(state)-1)];
    else                                          % Multiplicative Coupling
        neu_func = @(state, I) neu_model(state) * I;
    end
    
    
    I_osc = p.Results.ConstI;
    % Assume that the time is given in seconds, convert time to millisecond.
    dt = 1000*(t(2)-t(1));
    % Initialize membrane voltage vector.
    out = zeros(numel(t),numel(neu_state));
    % Simulate the neuron model.
    for i = 1:numel(I_ext)
        % Compute the derivatives of neuron state variables.
        d_neu_state = neu_func( neu_state, I_ext(i) );
        % Use the forward Euler method to integrate.
        neu_state   = neu_state + dt * d_neu_state;
        % Export the membrane voltage to output.
        out(i,:) = neu_state;
    end
    if ~p.Results.FullOutput,
        out = out(:,1);
    end


    % Define Default Hodgkin-Huxley Neuron ODEs. This function handle is 
    % written by Yiyin Zhou.
    function dx = hhn_ode(x)

        v = x(1);
        % compute dn
        a = exp(-(v+55)/10)-1;
        if a == 0
            dn = (1-x(2)) * 0.1 - x(2) * (0.125*exp(-(v+65)/80));
        else
            dn = (1-x(2)) * (-0.01*(v+55)/a) - x(2) * (0.125*exp(-(v+65)/80));
        end
        % compute dm
        a = exp(-(v+40)/10)-1;
        if a == 0
            dm = (1-x(3)) - x(3) * (4*exp(-(v+65)/18));
        else
            dm = (1-x(3)) * (-0.1*(v+40)/a) - x(3) * (4*exp(-(v+65)/18));
        end
        dh = (1-x(4)) * (0.07*exp(-(v+65)/20)) - x(4) / (exp(-(v+35)/10)+1);
        dV = I_osc - 120*x(3).^3*x(4)*(v-50) - 36 * x(2).^4 * (v+77) - 0.3 * (v+54.387);
        dx = [dV dn dm dh];
    end
    % Define Default Rinzel Neuron model ODEs.
    function dx = rinzel_ode(x)

        v = x(1);

        a = exp(-(v+55)/10)-1;
        alpha_n = ([a==0 a~=0])*[0.01; -0.01*(v+55)/a];
        a = exp(-(v+40)/10)-1;
        alpha_m = ([a==0 a~=0])*[1; -0.1*(v+40)/a];
        alpha_h = 0.07*exp(-(v+65)/20) ;                     

        beta_n = 0.125*exp(-(v+65)/80); 
        beta_m = 4*exp(-(v+65)/18);
        beta_h = 1/(exp(-(v+35)/10)+1);                  

        n_infty = alpha_n/(alpha_n + beta_n);
        m_infty = alpha_m/(alpha_m + beta_m);    
        h_infty = alpha_h/(alpha_h + beta_h);

        % compute dR
        S = 1.2714;  % s = (1-h_inf(0)) / n_inf(0) 
        R_infty = S/(1+S^2)*(n_infty + S*(1-h_infty));
        tau_R = 1 + 5*exp(-(v+55)^2/55^2);   
        dR = -3*(x(2)-R_infty) / tau_R;

        % Update the ionic currents and membrane voltage:
        dV = I_osc -120*m_infty^3*(1-x(2))*(v-50) - 36*(x(2)/S).^4*(v+77) - 0.3*(v+54.387);
        dx = [dV dR];
    end
    % Define Default Wilson Neuron model ODEs.
    function dx = wilson_ode(x)
        v = x(1);
        R_infty = 0.0135*v+1.03;
        dR = -(x(2)-R_infty) / 1.9;
        dV = I_osc -(17.81+0.4771*v+3.263e-3*v^2)*(v-55) - 26*x(2)*(v+92);  
        dx = [dV dR];
    end
end