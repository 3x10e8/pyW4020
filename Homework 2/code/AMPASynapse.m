%% AMPA Synapse taking Vpre as input
function [s,cond] = AMPASynapse(t,NT,params)

gmax = params.gmax;

dt = t(2)-t(1);
ddt = dt * 1e3; % change to msec

s = zeros(size(t));    % state
cond = zeros(size(t)); % conductance

for tt = 2:length(t)
    [s(tt),cond(tt)] = AMPASynapse_step(ddt,NT(tt-1),gmax,s(tt-1));
end

return