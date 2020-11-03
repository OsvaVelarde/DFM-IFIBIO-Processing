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
@title: xxx
"""

import numpy as np
import matplotlib.pyplot as plt

PATH = '/mnt/BTE2b/DBS/2020/GMdata/'

# Input file ------------------------------------------------
inputPATH = PATH + 'Biomk_method/Results/PSD/'
filename = inputPATH + 'PowerBands.dat'

pwbands = np.loadtxt(filename)

# Normalization of pwbands = 1 ------------------------------
#norm = pwbands[:,5:].sum(axis=1)
#pwbands[:,5:] = pwbands[:,5:] / norm[:, np.newaxis]
# ==========================================================

# Files Information 
indexfiles = ['001','002','003','004']
num_files = len(indexfiles)

# States 
state = ['Lesionado', 'Sham']
num_state = len(state)
color_state = ['bo','ro']

# Movement Information
mov_state =  ['Movement','Quiet']
num_movst = len(mov_state)
# ==========================================================

# Electrode Information
probe32 = np.loadtxt(PATH + 'Biomk_method/Files/Electrode_32.dat')
probe64 = np.loadtxt(PATH + 'Biomk_method/Files/Electrode_64.dat')
probes = [probe32,probe64,probe64,probe64]

# ==========================================================

# Frequency bands
lims_bands = [1,10,30,100,200,600]
num_bands = len(lims_bands)
name_band = ['0-1 Hz','1-10 Hz','10-30 Hz','30-100 Hz','100-200 Hz','200-600 Hz']

# ==========================================================
# Output parameters
outputPATH = PATH + 'Biomk_method/Results/PSD/Plots_PB/'

# ==========================================================
data_local = [0] * num_state

for ff in range(num_files):

	# Electrode - Structure -----------------------------------
	dimy, dimx = np.shape(probes[ff])

	if dimx >= dimy:
		dimh = dimx
		dimv = dimy
		order='C'
	else:
		dimh = dimy
		dimv = dimx
		order='F'

	# Order of data for plotting -----------------------------
	permutation = (np.reshape(probes[ff],(dimx*dimy,),order=order)-1).astype(int)

	for mm in range(num_movst):

		# Plot I: Map between 2 bands
		FIGI, AXSI = plt.subplots(1,2)

		# Selection of data: Sham vs Lesionado ------------------
		for ss in range(num_state):
			condition = (pwbands[:,0] == ff) & (pwbands[:,1] == ss) & (pwbands[:,2] == mm) 
			data_local[ss] = pwbands[condition,4:]
			data_local[ss] = data_local[ss][permutation,:]

			# Plot I ---------------------------------------------
			for bb in range(2):
				AXSI[bb].plot(data_local[ss][:,1],data_local[ss][:,2+bb],color_state[ss])
				AXSI[bb].set_xlabel(name_band[1])
				AXSI[bb].set_ylabel(name_band[2+bb])
				AXSI[bb].set_xlim([0,0.1])
		
		figname = outputPATH+ 'Maps_File_' + indexfiles[ff] + '_' + mov_state[mm] + '.eps'
		FIGI.savefig(figname,format='eps')
		plt.close(FIGI)

		# Plot II: For (file, mov_state, band) fixed ---------------
		for bb in range(num_bands):

			FIGII, AXSII = plt.subplots(dimv,1)
			FIGII.suptitle(name_band[bb])
			AXSII[dimv-1].set_xlabel('Canal')

			for vv in range(dimv):
				for ss in range(num_state):
					AXSII[vv].plot(data_local[ss][dimh*vv:dimh*(vv+1),bb],color_state[ss])
					AXSII[vv].set_xticks(range(dimh))

			figname = outputPATH+ 'File_' + indexfiles[ff] + '_' + mov_state[mm] + '_' + name_band[bb] + '.eps'
			FIGII.savefig(figname,format='eps')
			plt.close(FIGII)