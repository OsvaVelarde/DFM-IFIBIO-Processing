# -*- coding: utf-8 -*-
"""
Biomarkers in Parkinson's Disease (PD): Cross Frequency Coupling.

Description:
Analysis of the data measured in Gustavo Murer's laboratory.

Data: LFP measured in mice in parkinsonian and control states.
Reference:

Analysis: PSD calculation and spectrogram.

@author: Osvaldo M Velarde
@date: 30 June 2020
@title: Biomarkers in PD.
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

import sklearn.preprocessing
from sklearn.preprocessing import scale

from scipy.signal import iirnotch
from scipy.signal import butter
from scipy.signal import filtfilt
from scipy.signal import periodogram
from scipy.signal import freqz

sys.path.append('/mnt/BTE2b/DBS/2020/GMdata/Biomk_method/Scripts/Modules/')

import comodulogram
import filtering
from segmentation import function_segmentation
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

chnshow = [[16], range(64), range(64), range(64)]
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

# Windows parameters
name_window='hann'
sflag = True # %True = 'symmetric' - False = 'periodic'
window_parameters = {'name':name_window, 'sflag':sflag}   

fs = 1250 # [Hz]
pad = 0
method_parameters = {'Pad':pad,'fs':fs}

print('Parámetros cargados')

# ==========================================================

# Outputs parameters ---------------------------------------
#outputPATH = PATH + 'Biomk_method/Results/PSD/'
#output_filename = outputPATH + 'PowerBands.dat'
#output_file = open(output_filename,'a')

# ==========================================================

# Pre-processing parameters.
freq_notch = 50.3 
Q_notch = 30

bNotch,aNotch = iirnotch(2*freq_notch/fs,Q_notch)

# ==========================================================

# Bands of interest: Frequency bands of interest for  PLV computation. # [Hz].
freq_bands = [[0.05, 3], [30, 150]]            
name_bands = {'1-19 Hz','20-200 Hz'};

fXres = 1
fYres = 1

fXmin = (freq_bands[0][1]+freq_bands[0][0])/2
fXmax = fXmin + fXres / 2
fXBw  = freq_bands[0][1]-freq_bands[0][0]

fYmin = (freq_bands[1][1]+freq_bands[1][0])/2
fYmax = fYmin + fYres / 2 
fYBw  = freq_bands[1][1]-freq_bands[1][0]

# Comodulogram. 
# For the description of the parameters see function "function_CFCcfg".            

BPFXcfg = {'Bw':fXBw,
           'zeropadding':0,
           'freqWindowParam': window_parameters,
           'timeWindowParam': window_parameters,
           'conv':'circular',
           'causal':0,
           'Nf':1,
           'function':'function_FDF'}

BPFYcfg = {'Bw':fYBw,
           'zeropadding':0,
           'freqWindowParam': window_parameters,
           'timeWindowParam': window_parameters,
           'conv':'circular',
           'causal':0,
           'Nf':1,
           'function':'function_FDF'}

CFCcfg = {'fXmin':fXmin, 'fXmax':fXmax, 'fXres':fXres,
          'fYmin':fYmin, 'fYmax':fYmax, 'fYres':fYres,
          'fXlookAt':'PHASE', 'fYlookAt':'PHASEofAMPLITUDE',
          'nX':1, 'nY':1,
          'BPFXcfg':BPFXcfg, 'BPFYcfg':BPFYcfg,
          'saveBPFsignal':0,                 
          'Nbins':18, 'sameNumberOfCycles':0,
          'CFCmethod':'plv', 'verbose':1,
          'perMethod':'FFTphaseShuffling',
          'Nper':100, 'Nrep':1, 'Pvalue':0.05,
          'corrMultComp':'Bonferroni',
          'fs':fs}

CFCcfg = comodulogram.function_setCFCcfg(CFCcfg) 

print(CFCcfg['fXcfg']['BPFcfg'])
exit()
# %% Movement vs. Quiet information
# up_bound = ones(1,100);             

# %% Checking filter flags.
CHEKING_FILTER_FLAG = 1;
plotFlag = 0;
            
# %% Input parameters to compute the comodulogram.
segLen = 200 # [sec]
peroverlap = 90 # %
        
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
    #index_artifacts, signal = preproc.remove_artifacts(signal,lenwinRA,coeffRA)

    # Pre-processing: Eliminar 50 Hz?
    signal = filtfilt(bNotch,aNotch,signal,axis=0)

    # Pre-processing: Referecing
    signal = preproc.referencing(signal,ref_cfg)

    # Selection of channels for show. 
    signal = signal[:,chnshow[ff]]

    # Pre-processing: Normalization
    signal = sklearn.preprocessing.scale(signal)
    Nsamples = np.shape(signal)[0] # Compute the number of samples.

    print('Datos pre-procesados')

    # =============================================

    # # Segmentación de periodos: Quieto - Movimiento        
    # indexmov = [[],[]]

    # for mm in range(num_movst):
    #   filename =  prefname + '-IndicesLFP_' + mov_state[mm] + 'Periods.dat'
    #   indexmov[mm] = np.loadtxt(filename,dtype='int')

    # print('Datos sobre Movimiento cargados')
    
    # =============================================
    # Check the filters and compute the settling time (percLevel).

    indSettling = 0 # Initialize the index.             

    #  Check the filters just once.
    if CHEKING_FILTER_FLAG == 1:
      CHEKING_FILTER_FLAG = 0

      # BPFs for the "x-y" axis of the comodulogram.
      indSettlingLF = filtering.function_checkFilter(CFCcfg['fXcfg']['BPFcfg'], CFCcfg['fs'], Nsamples, plotFlag);
      indSettlingHF = filtering.function_checkFilter(CFCcfg['fYcfg']['BPFcfg'], CFCcfg['fs'], Nsamples, plotFlag);

      # Compute the maximum settling time.
      indSettling = np.amax(np.concatenate((indSettlingHF,indSettlingLF.T)))

    # Compute the maximum settling time.
    indSettling = max(indSettling, Nsamples+1)
 
    print('Filtros chequeados')

    # =============================================
    # Loops across the channels.

    for kk in range(num_chnshow):
      print('Canal ',kk)

      signal_channel = signal[:,kk]

      # =============================================
      # Reflect the time series to minimize edge artifacts ----------------- 
      # due to the transient response of the BPFs. -------------------------

      signal_channel = np.concatenate((signal_channel[::-1],signal_channel,signal_channel[::-1]))
      signal_channel = signal_channel.reshape(signal_channel.shape[0],-1)    

      # Compute the Band-Pass Filtering. -----------------------------------

      LFSIGNAL, _ = comodulogram.function_comodulogramBPF(signal_channel, CFCcfg['fXcfg']['BPFcfg'], CFCcfg['fs'], indSettling)
      HFSIGNAL, _ = comodulogram.function_comodulogramBPF(signal_channel, CFCcfg['fYcfg']['BPFcfg'], CFCcfg['fs'], indSettling)

      # Compute the length of "LFSIGNAL" to compute the number of segments.
      N_LFSIGNAL = np.shape(LFSIGNAL)[0]
      N_HFSIGNAL = np.shape(HFSIGNAL)[0]

      # Restore the length of the raw signal.
      #signal_channel = signal_channel[indSettling-1:signal_channel.shape[0]-(indSettling-1),:]

      # =============================================
      # Phase - PhaseofAmplitude of LF and HF (all signal)

      # # LF segment ----------------             
      # for jj in range(LFSIGNAL.shape[2]):
      #   LFSIGNAL[:,:,jj] = scale(LFSIGNAL[:,:,jj],axis=0) #%Z-score normalization of the segment

      # # Reflect the segment to minimize edge artifacts due to the transient response of the Hilbert transform.
      # LFSIGNAL = np.concatenate((LFSIGNAL[::-1,:,:],LFSIGNAL,LFSIGNAL[::-1,:,:]))

      # # Compute the phase and phase of amplitude
      # PhaseLF, _ = comodulogram.function_comodulogramFeature(LFSIGNAL, CFCcfg['fXcfg'], CFCcfg['fs'], N_LFSIGNAL+1)
      # PhaseLF = np.squeeze(PhaseLF,axis=2)

      #fXcfg_aux = CFCcfg['fXcfg']
      #fXcfg_aux['lookAt'] = 'PHASEofAMPLITUDE'
      # PhaseofAmplitudeLF, _ = comodulogram.function_comodulogramFeature(LFSIGNAL, fXcfg_aux, CFCcfg['fs'], N_LFSIGNAL+1)
      # PhaseofAmplitudeLF = np.squeeze(PhaseofAmplitudeLF,axis=2)
      
      # # HF segment -----------------
      # for jj in range(HFSIGNAL.shape[2]): #%Loop for Bandwidths.
      #   HFSIGNAL[:,:,jj] = scale(HFSIGNAL[:,:,jj],axis=0)

      # # Reflect the segment to minimize edge artifacts due to the transient response of the Hilbert transform.
      # HFSIGNAL = np.concatenate((HFSIGNAL[::-1,:,:],HFSIGNAL,HFSIGNAL[::-1,:,:]))

      # # Compute the phase/amplitude/frequency segment.
      # PhaseofAmplitudeHF, _ = comodulogram.function_comodulogramFeature(HFSIGNAL, CFCcfg['fYcfg'], CFCcfg['fs'], N_HFSIGNAL+1)
      # PhaseofAmplitudeHF = np.squeeze(PhaseofAmplitudeHF,axis=2)

      #fYcfg_aux = CFCcfg['fYcfg']
      #fYcfg_aux['lookAt'] = 'PHASE'
      # PhaseHF, _ = comodulogram.function_comodulogramFeature(HFSIGNAL, fYcfg_aux, CFCcfg['fs'], N_HFSIGNAL+1)
      # PhaseHF = np.squeeze(PhaseHF,axis=2)

      # # =============================================

      # print(np.shape(PhaseLF),np.shape(PhaseHF),np.shape(PhaseofAmplitudeLF),np.shape(PhaseofAmplitudeHF))

      # FIG, AXS = plt.subplots()
      # AXS.plot(PhaseLF)
      # AXS.plot(PhaseofAmplitudeLF)
      # AXS.plot(PhaseHF)
      # AXS.plot(PhaseofAmplitudeHF)

      # =============================================

      # Segmentation of the Band-Pass Filtered time series.
      indSegment = function_segmentation(segLen, peroverlap, N_LFSIGNAL, CFCcfg['fs'])
      metrics = np.zeros((len(indSegment),3))

      print('Señal fragmentada')

      for ii in range(len(indSegment)): # Loop across the segments.

        print('Segmento ',ii,'de ',len(indSegment))
        # Compute the LF/HF segment ---------------------------------------------------

        segLFSIGNAL = np.copy(LFSIGNAL[indSegment[ii,0]:indSegment[ii,1],:,:])
        segHFSIGNAL = np.copy(HFSIGNAL[indSegment[ii,0]:indSegment[ii,1],:,:])
        Namples_per_segment = len(segLFSIGNAL) # Compute the number of samples of the segment.

        # ============================================= 

        # Compute the phase/amplitude/frequency for the LF segment ----------------             
        for jj in range(segLFSIGNAL.shape[2]):
          segLFSIGNAL[:,:,jj] = scale(segLFSIGNAL[:,:,jj],axis=0) #%Z-score normalization of the segment

        # Reflect the segment to minimize edge artifacts due to the transient response of the Hilbert transform.
        segLFSIGNAL = np.concatenate((segLFSIGNAL[::-1,:,:],segLFSIGNAL,segLFSIGNAL[::-1,:,:]))

        # Compute the phase/amplitude/frequency segment.
        segPLF, _ = comodulogram.function_comodulogramFeature(segLFSIGNAL, CFCcfg['fXcfg'], CFCcfg['fs'], Namples_per_segment+1)
        segPLF = np.squeeze(segPLF,axis=2)

        fXcfg_aux = CFCcfg['fXcfg'].copy()
        fXcfg_aux['lookAt'] = 'PHASEofAMPLITUDE'

        segPALF, _ = comodulogram.function_comodulogramFeature(segLFSIGNAL, fXcfg_aux, CFCcfg['fs'], Namples_per_segment+1)
        segPALF = np.squeeze(segPALF,axis=2)

        # --------------------------------------------------------------------------

        # Compute the phase/amplitude/frequency for the HF segment -----------------
        for jj in range(segHFSIGNAL.shape[2]): #%Loop for Bandwidths.
          segHFSIGNAL[:,:,jj] = scale(segHFSIGNAL[:,:,jj],axis=0)

        # Reflect the segment to minimize edge artifacts due to the transient response of the Hilbert transform.
        segHFSIGNAL = np.concatenate((segHFSIGNAL[::-1,:,:],segHFSIGNAL,segHFSIGNAL[::-1,:,:]))

        # Compute the phase/amplitude/frequency segment.
        segPAHF, _ = comodulogram.function_comodulogramFeature(segHFSIGNAL, CFCcfg['fYcfg'], CFCcfg['fs'], Namples_per_segment+1)

        fYcfg_aux = CFCcfg['fYcfg'].copy()
        fYcfg_aux['lookAt'] = 'PHASE'

        segPHF, _ = comodulogram.function_comodulogramFeature(segHFSIGNAL, fYcfg_aux, CFCcfg['fs'], Namples_per_segment+1)

        # ============================================= 

        # Compute the Phase Locking Value ------------------------------------------
        PPC, _, _ = comodulogram.function_PLV(segPLF, segPHF, 0, 0, CFCcfg)
        PAC, _, _ = comodulogram.function_PLV(segPLF, segPAHF, 0, 0, CFCcfg)
        AAC, _, _ = comodulogram.function_PLV(segPALF, segPAHF, 0, 0, CFCcfg)

        metrics[ii,0] = np.abs(PPC)
        metrics[ii,1] = np.abs(PAC)
        metrics[ii,2] = np.abs(AAC)


      FIG2, AXS2 = plt.subplots(1,3)
      AXS2[0].plot(metrics[:,0])
      AXS2[1].plot(metrics[:,1])
      AXS2[2].plot(metrics[:,2])

      plt.show()
      exit()      
