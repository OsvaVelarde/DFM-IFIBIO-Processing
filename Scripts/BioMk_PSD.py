# -*- coding: utf-8 -*-
"""
Biomarkers in Parkinson's Disease (PD): Espectral analysis.

Description:
Analysis of the data measured in Gustavo Murer's laboratory.

Data: LFP measured in mice in parkinsonian and control states.
Reference:

Analysis: PSD calculation and spectrogram.

@author: Osvaldo M Velarde
@date: 13 April 2020
@title: Dynamical Bayesian Inference in Van der Pol Oscillator
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

import sklearn.preprocessing

from scipy.signal import get_window
from scipy.signal import periodogram
from scipy.signal import spectrogram

sys.path.append('/mnt/BTE2b/DBS/2020/GMdata/Biomk_method/Scripts/Modules/')

from biomarkers import nextpow2
from biomarkers import power_in_bands
import preprocessing as preproc

PATH = '/mnt/BTE2b/DBS/2020/GMdata/'

# Data sets ------------------------------------------------
dataPATH = PATH + 'data/'

# States & Files Information
state = ['Lesionado/20170418-', 'Sham/20170216-']
num_state = len(state)

#indexfiles = ['001']
indexfiles = ['001','002','003','004']
num_files = len(indexfiles)

# Channels Information
num_channels = [32,64,64,64]
MAXChn = max(num_channels)

# Electrode Information
probe32 = np.loadtxt(PATH + 'Biomk_method/Files/Electrode_32.dat')
probe64 = np.loadtxt(PATH + 'Biomk_method/Files/Electrode_64.dat')
probes = [probe32,probe64,probe64,probe64]

chnshow = [range(32), range(64), range(64), range(64)]
# ShowChannels = {17, 36, 35:38, 36};
# %ShowChannels = {[7 23 5 21 3 19 1 17], [4,5,12,13,20,21,28,29,36,37,44,45,52,53,60,61], [1 8 9 16 17 24 25 32 33 40 41 48 49 56 57 64], [1 8 9 16 17 24 25 32 33 40 41 48 49 56 57 64]};
# %ShowChannels = {reshape(Probe32',[1,32]), reshape(Probe64',[1,64]), reshape(Probe64',[1,64]), reshape(Probe64',[1,64])};

# Movement Information
mov_state =  ['Movement','Quiet']
num_movst = len(mov_state)
mintimeMQ = 3 # [sec]

# ==========================================================

# Parameters of Processing ---------------------------------

ref_cfg = {'method':'hard','electrode':0}

# Frequency bands
lims_bands = [1,10,30,100,200,600]
num_bands = len(lims_bands)

# Windows parameters
name_window='hann'
sflag = True # %True = 'symmetric' - False = 'periodic'
window_parameters = {'name':name_window, 'Alpha':sflag}   

fs = 1250 # [Hz]
pad = 0
method_parameters = {'Pad':pad,'fs':fs}

# Method of spectrum.
conflevel = 0.95
onesided = True

print('Parámetros cargados')

# ==========================================================

# Outputs parameters ---------------------------------------
outputPATH = PATH + 'Biomk_method/Results/PSD/'
output_filename = outputPATH + 'PowerBands.dat'
output_file = open(output_filename,'a')

# ==========================================================

for ff in range(num_files):

	num_chnshow = len(chnshow[ff])
	ref_cfg['electrode'] = probes[ff]

	for ss in range(num_state):
		
		prefname = dataPATH + state[ss] + indexfiles[ff]
		
		# Load data	
		filename = prefname + '.dat'
		signal = np.loadtxt(filename)[:,0:num_channels[ff]]

		print('Datos leidos')

		# =============================================

		# Pre-processing: Remove artifacts.
		lenwinRA = 6250
		coeffRA = 0.62 # Cambiar si es necesario
		index_artifacts, signal = preproc.remove_artifacts(signal,lenwinRA,coeffRA)

		# Pre-processing: Referecing
		signal = preproc.referencing(signal,ref_cfg)

		# Selection of channels for show. 
		signal = signal[:,chnshow[ff]]

		# Pre-processing: Normalization
		signal = sklearn.preprocessing.scale(signal)


		print('Datos pre-procesados')
				
		# =============================================

		# Segmentación de periodos: Quieto - Movimiento        
		indexmov = [[],[]]
		frags_signals = [[],[]] #frags_signal[mov][periods][time,channels]

		for mm in range(num_movst):
			filename =  prefname + '-IndicesLFP_' + mov_state[mm] + 'Periods.dat'
			indexmov[mm] = np.loadtxt(filename,dtype='int')
			element, _, _ = preproc.fragmentation(signal,indexmov[mm],mintimeMQ*fs)
			frags_signals[mm] = element

		num_fragms = [len(frags_signals[mm]) for mm in range(num_movst)]

		del signal
		print('Datos fragmentados')
		# ===================================================

		for mm in range(num_movst):

			PSD = [0] * num_chnshow			# Array-2d: (num_channels x num_freqs)

			# For particular (file, state, movstate) computing the PSD for all channels and fragments. 
			aux_fragms = [nextpow2(np.shape(frags_signals[mm][jj])[0]) for jj in range(num_fragms[mm])]
			maxaux = max(aux_fragms)

			for jj in range(num_fragms[mm]):
				num_samples = np.shape(frags_signals[mm][jj])[0]

				pad = maxaux - aux_fragms[jj]
				nfft = 2^(nextpow2(num_samples) + pad)
				window = get_window(window_parameters['name'], num_samples)

				for kk in range(num_chnshow):
					#fPSD, PSD[jj][kk] = periodogram(frags_signals[mm][jj][:,kk],method_parameters['fs'],window,nfft,return_onesided=onesided)
					fPSD, PSD_jj_kk = periodogram(frags_signals[mm][jj][:,kk],method_parameters['fs'],window,nfft,return_onesided=onesided)
					print(np.shape(PSD_jj_kk))
					PSD[kk] = PSD[kk] + PSD_jj_kk

			PSD = np.array(PSD)/num_fragms[mm] 			# Average over number of fragments

			# ===================================================

			# Processing of PSD in each channels
			for kk in range(num_chnshow):
				
				# Power in frequency bands 
				pwbands = power_in_bands(fPSD,PSD[kk],lims_bands)
				cte_normalization = np.sum(pwbands)

				# Save data
				output = np.reshape(np.concatenate(([ff,ss,mm,kk],pwbands)),(1,4 + num_bands))
				np.savetxt(output_file,output)

				# Normalization of PSD using Integral
				PSD[kk] = PSD[kk]/cte_normalization 
					
				# Normalization of PSD using median (Implementar)
				# %median_psd = ones(size(Mean_PSD_Quiet,1),1) * median(Mean_PSD_Quiet);
				# %Mean_PSD_Quiet = Mean_PSD_Quiet ./ median_psd;
                             
			# ===================================================

			# Average over number of channels
			meanPSD = np.mean(PSD,axis=0) #Array-1d: (num_freqs)

			# ===================================================

			# Plots PSDs
			FIG,AXS = plt.subplots(1,1)
			AXS.plot(fPSD,10*np.log10(PSD.T),'b',alpha=0.5)
			AXS.plot(fPSD,10*np.log10(meanPSD),'r')

			AXS.set_xlim([0.01,650])
			AXS.set_ylim([-80,40])
			AXS.set_xlabel('Frecuencia (Hz)')
			AXS.set_ylabel('PSD (dB)')
			AXS.set_xscale('log')

			figname = outputPATH + 'Plots_PSD/' + 'PSD_State_' + str(ss) + '_File_'+ indexfiles[ff] + '_' + mov_state[mm] + '.eps'
			FIG.savefig(figname,format='eps')			

output_file.close()