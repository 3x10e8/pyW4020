%% NMVDA Synapse 
% NMDA Synapse is affected by Mg2+ block which is dependent on voltage of
% postsynaptic neuron
function [s,cond] = NMDASynapse(t,NT,Vpost,Mg,params)

gmax = params.gmax;

dt = t(2)-t(1);
ddt = dt * 1e3; % change to msec

s = zeros(size(t));    % state
cond = zeros(size(t)); % conductance

for tt = 2:length(t)
    [s(tt),cond(tt)] = NMDASynapse_step(...
        ddt,NT(tt-1),Vpost(tt-1),Mg(tt-1),gmax,s(tt-1));
end

return