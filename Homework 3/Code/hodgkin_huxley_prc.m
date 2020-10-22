%HODGKIN_HUXLEY_PRC find the phase response curve of a given neuron model
%   [PERIOD LIMITCYCLE PRC] = HODGKIN_HUXLEY_PRC(DT,I_OSC,VARARGIN) finds
%   the phase response curve PRC of a given neuron model with time step DT
%   and constant current bias I_OSC. The limit cycle, LIMITCYCLE, and its
%   period PERIOD are computed beforehand in order to find the PRC. Note
%   that PERIOD is the length of the limit cycle.
%
%   [...] = HODGKIN_HUXLEY_PRC(...,'Method',METHOD) finds PRC using 
%   specified METHOD. When METHOD is set to 'Malkin' as default, the 
%   Malkin's method is applied. When 'Impulse' is specified, the given
%   neuron model is trigger with current impulse for each state of limit
%   cycle, and the resultant phase shift is measured.
%   
%   [...] = HODGKIN_HUXLEY_PRC(...,'Model',MODEL) finds PRC for the given
%   neuron model MODEL. MODEL can be one of 'HH', as default, 'Rinzel', and
%   'Wilson', or a Matlab function handle for user defined model. When user 
%   defined model is provided, please also use 
%   HODGKIN_HUXLEY_PRC(...,'State',state) to specify the initial values of
%   the state variables for the given model.
%
%   Author:      Chung-Heng Yeh <chyeh@ee.columbia.edu>
%
%   Copyright 2010-2012   Chung-Heng Yeh

function [period limitCycle PRC] = hodgkin_huxley_prc(dt,I_osc,varargin)    
    % Handle the optional input parameters.
    % =====================================================================
    p = inputParser;
    p.KeepUnmatched = true;
    % Specify the Neuron ODEs.
    addParamValue(p,'Model','HH',@(x) isa(x,'function_handle') || any(validatestring(x,{'HH','Wilson','Rinzel'})));
    % Specify the initial values of neuron state variables.
    addParamValue(p,'State',0, @isnumeric);
    % Specify the method for computing PRC.
    addParamValue(p,'Method','Malkin',@(x) any(validatestring(x,{'Malkin','Winfree'})));

    % Parse the optional parameters.
    parse(p,varargin{:});

    % Handle the unexpected input parameter.
    UnmatchedParam = fieldnames(p.Unmatched);
    if ~isempty(UnmatchedParam)
        error(['"',UnmatchedParam{1},'" is not a valid parameter.']);
    end

    % Define spike detection function handle.
    spike_detect = @(V) (V(1) < V(2) & V(3) < V(2));
    
    % Setup the numerical method.
    % =====================================================================
    % Initialize the Neuron function handel and state variables. 
    if strcmpi( p.Results.Model, 'HH' ),           % Hodgkin-Huxley Model
        neu_func = @hhn_ode;
        neu_state = [0 0 0 1];
    elseif strcmpi( p.Results.Model, 'Rinzel' ),   % Rinzel Model
        neu_func = @rinzel_ode;
        neu_state = [-60 0.4];        
    elseif strcmpi( p.Results.Model, 'Wilson' ),   % Wilson Model
        neu_func = @wilson_ode;
        neu_state = [-70 -0.088];        
    else                                           % User specified Model
        neu_func = p.Results.Model;
        neu_state = p.Results.State;         
    end
    
    % Assume that the time is given in seconds, convert time to millisecond.
    dt = 1000*dt;
    
    % Find the limit cycle.
    % =====================================================================
    neu_state_at_spike = Inf(size(neu_state));
    prev_neu_state = neu_state;
    counter = 0;
    spikes = 0;
    while spikes<8
        next_neu_state = neu_state + dt * neu_func( neu_state );
        if spike_detect( [prev_neu_state(1) neu_state(1) next_neu_state(1)] )
            neu_state_at_spike = neu_state;
            period  = counter;
            counter = 0;
            spikes = spikes + 1;
            %disp(spikes);
        end
        % Update the neuron state variables.
        prev_neu_state = neu_state;
        neu_state = next_neu_state;
        counter = counter + 1;
        if counter > 2e7
            error('Limit Cycle does not exist. Use smaller time step...');
        end
    end
    % Record the limit cycle.
    counter = period;
    neu_state = neu_state_at_spike;
    limitCycle = zeros(period,numel(neu_state));
    for i = 1:counter
        limitCycle(i,:) = neu_state;
        neu_state = neu_state + dt * neu_func( neu_state );
    end
    
    % Find the phase response curve.
    % =================
    % Initialize the PRC matrix.
    PRC = zeros(period,numel(neu_state));
    
    if strcmpi(p.Results.Method,'Malkin')
        % Malkin's Method - This method is originally written by Yiyin Zhou
        
        % Initialize the Jacobian of the model
        A  = zeros( numel(neu_state) );
        % Initialize the perturbation
        dx = 1e-5;
        dX = dx*eye(numel(neu_state)); 

        PRC(1,:) = neu_func(neu_state);
        for iters=1:50
            PRC(end,:) = PRC(1,:);
            for i = period:-1:2
                for j = 1:numel(neu_state)
                    A(j,:) = (neu_func(limitCycle(i,:)+dX(j,:)) - neu_func(limitCycle(i,:)))/dx;
                end
                PRC(i-1,:) = PRC(i,:) + dt*(PRC(i,:)*A');
            end
        end
        PRC = PRC/(PRC(1,:)*neu_func(limitCycle(1,:))');
    else
        % Winfree's Method - trigger the neuron with current impulse and 
        % measure the phase shift
        spike_period = period*dt;
        step = 50;
        for i = 1:step:period
            t = (i-1)*dt;
            neu_state = limitCycle(i,:);
            neu_state(1,1) = neu_state(1,1) + 1;
            prev_neu_state = neu_state;
            spike_count = 0;
            while spike_count < 3,
                t = t+dt;
                next_neu_state = neu_state + dt * neu_func( neu_state );
                if spike_detect( [prev_neu_state(1) neu_state(1) next_neu_state(1)] )
                    spike_count = spike_count + 1;
                end
                % Update the neuron state variables.
                prev_neu_state = neu_state;
                neu_state = next_neu_state;
            end
            PRC(i:i+step) = period*3*dt-t;
            if PRC(i) > spike_period/2, PRC(i:i+step) = PRC(i:i+step) - spike_period; end
        end
    end
   
    % Define Default Hodgkin-Huxley Neuron ODEs.
    function dx = hhn_ode(x)  
        % This function handle is written by Yiyin Zhou.
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
    % Define Default Rinzel Neuron ODEs.
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
    % Define Default Wilson Neuron ODEs.
    function dx = wilson_ode(x)
        v = x(1);
        R_infty = 0.0135*v+1.03;
        dR = -(x(2)-R_infty) / 1.9;
        dV = I_osc -(17.81+0.4771*v+3.263e-3*v^2)*(v-55) - 26*x(2)*(v+92);  
        dx = [dV dR];
    end
end
