classdef utils
    methods(Static)
        function [spike_states] = poissonSpikes(t, rate, nTrials)
            if nargin < 3
                nTrials = 1;
            end
            dt = diff(t(1:2)); % unit in s
            spike_states = rand([nTrials, length(t)]) < rate*dt;
        end
       
        function plotRaster(t, spike_state,cs)
            spike_state = logical(spike_state);
            nTrials = size(spike_state,1);
            if nargin<3
                cs = zeros(nTrials,3);
            end
            hold all;
            for tr = 1:nTrials
                tk = t(spike_state(tr,:)); tk = tk(:);
                plot([tk,tk]',[ones(size(tk))*tr-0.4,ones(size(tk))*tr+0.4]','k');
                plot([min(t),min(t)], [tr-0.5,tr+0.5],'Color',cs(tr,:),'LineWidth',3);
            end
            ylim([0,size(spike_state, 1)+1]);
            xlim([min(t),max(t)]);
        end
    end
end