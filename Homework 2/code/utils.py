import numpy as np
import matplotlib.pyplot as plt 

class utils:
	def poissonSpikes(t, rate, nTrials = 1):
		dt = t[1]-t[0] # unit in s
		spike_states = np.random.rand(nTrials, len(t)) < rate*dt
		return spike_states

	def plotRaster(t, spike_state, cs=[]):
        #spike_state = logical(spike_state)

		nTrials = spike_state.shape[0] # number of rows
		#print(nTrials)

		if cs == []: # define cs if not provided by user
			cs = np.zeros(shape=(nTrials,3))

		for tr in range(nTrials):
			tk = t[spike_state[tr,:]] # spike times
			print(tk)
			#tk = tk(:)

			#plt.plot([tk,tk]',[ones(size(tk))*tr-0.4,ones(size(tk))*tr+0.4]','k')
			plt.plot(tk, np.ones_like(tk)*tr -0.4, 'k')
			plt.plot(tk, np.ones_like(tk)*tr +0.4, 'k')
			#plt.plot([min(t),min(t)], [tr-0.5,tr+0.5],'Color',cs(tr,:),'LineWidth',3)

		'''
		plt.ylim([0,size(spike_state, 1)+1])
		plt.xlim([min(t),max(t)])
		'''
		plt.show()

if __name__ == '__main__': # run self-test

	t = np.arange(0, 1, 1e-3) # 1s time array in 1ms steps
	rate = 100 # seems to be Hz
	nTrials = 5
	spike_states = utils.poissonSpikes(t, rate, nTrials)

	# print spike train
	print(spike_states)

	# number of spikes by row, where each row is a trial
	print('#spikes: {} in {} seconds'.format(np.sum(spike_states, axis=1), t[-1]))

	utils.plotRaster(t, spike_states)