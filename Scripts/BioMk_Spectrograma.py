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
import preprocessing as preproc

# % Windows VERSION (Osva Notebook)
# %pathFunctions3 = 'C:\Users\osva_\Documents\Toolboxs_MATLAB\Espectrograma\';
# %pathFunctions5 = 'C:\Users\osva_\Documents\Toolboxs_MATLAB\Spectrum\';
# %pathFunctions6 = 'C:\Users\osva_\Documents\Toolboxs_MATLAB\Filtering\';

# % Windows VERSION (Osva Notebook)
# %pathFunctions1 = 'C:\Users\osva_\Documents\Toolboxs_MATLAB\buzcode-master\externalPackages\FMAToolbox\IO\';
# %pathFunctions2 = 'C:\Users\osva_\Documents\Toolboxs_MATLAB\buzcode-master\externalPackages\FMAToolbox\Helpers\';
# %pathFunctions3 = 'C:\Users\osva_\Documents\Toolboxs_MATLAB\Espectrograma\';
# %pathFunctions4 = 'C:\Users\osva_\Documents\Toolboxs_MATLAB\Normalization\';
# %pathFunctions5 = 'C:\Users\osva_\Documents\Toolboxs_MATLAB\Spectrum\';
# %pathFunctions6 = 'C:\Users\osva_\Documents\Toolboxs_MATLAB\Filtering\';


PATH = '/mnt/BTE2b/DBS/2020/GMdata/'

# Data sets ------------------------------------------------
dataPATH = PATH + 'data/'

# States & Files Information
state = ['Lesionado/20170418-', 'Sham/20170216-']
num_state = len(state)

indexfiles = ['001','002','003','004']
num_files = len(indexfiles)

# Channels Information
num_channels = [32,64,64,64]
MAXChn = max(num_channels)

# Electrode Information
probe32 = np.loadtxt(PATH + 'Biomk_method/Files/Electrode_32.dat')
probe64 = np.loadtxt(PATH + 'Biomk_method/Files/Electrode_64.dat')
probes = [probe32,probe64,probe64,probe64]

chnshow = [[16], [36], range(35,38), [36]]

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
#color_band = ['red','blue','black','cyan','green']
legend_band = ['0-1 Hz','1-10 Hz','10-30 Hz','30-100 Hz','100-200 Hz','200-600 Hz']


# Pre-filtering parameters
order = 4
timesf = 1
fc =  0.8
zeroPhase = 1
percLevel = 1
hpf_cfg = {'order':order, 'times':timesf , 'fc':fc,
			'zeroPhase':zeroPhase, 'percLevel':percLevel,
			'function':'function_butterHPF_v1'}


# Method of spectrogram.
#correctionAmp = 'z-score'
# NW = 3;
# DropLastTaper = 0; 
# MTMmethod = 'eigen'; % Other option: 'adapt'
name_window='hann'
pad = 0
fs = 1250 					# [Hz]
segmentwitdh = 5 			# [sec]
overlapper = 0.9 			# (0,1)
nperseg = segmentwitdh * fs # int
noverlap = int(0.9*nperseg)
onesided = True

# Filters-Hilbert method

# BPFcfg = struct('f1',LimBands(:,1),...
#                 'f2',LimBands(:,2),...
#                 'freqWindowParam',windowParam,...
#                 'conv','circular',... 
#                 'causal',0,...
#                 'zeropadding',0,...
#                 'timeWindowParam',windowParam,...
#                 'function','function_FDF');            
            
# SFHcfg = struct('BPFcfg',BPFcfg,...
#                 'correctionAmp',correctionAmp,...
#                 'freqMarkers',[],... 
#                 'textMarkers',[],... 
#                 't0', 0,...
#                 'fs',fs,...
#                 'plot',0);   
   
# Parameters of Plotting ----------------------------------


# ==========================================================

print('Parámetros cargados')

# %% Pre-allocation
# PowerBandsQuiet = zeros(MAXChn,NBand,NumberStates,NNf);
# pwband_quiet = np.zeros()
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

		# Pre-processing: Espectrograma.
		#[HPFsignal, HPFmag, f_aux] = function_butterHPF_v1(signal, HPFcfg, fs);
		#indSettling = size(T,1); # Filter Hilbert Parameters 	
		#signal = HPFsignal;
		#clear HPFsignal HPFmag

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

		for mm in range(num_movst):
			filename =  prefname + '-IndicesLFP_' + mov_state[mm] + 'Periods.dat'
			indexmov[mm] = np.loadtxt(filename,dtype='int')

		print('Datos sobre Movimiento cargados')
		# =============================================

		# Espectrograma ---------------------------------

		num_samples = np.shape(signal)[0]
		nfft = 2^(nextpow2(num_samples) + pad)

		f, t, Sxx = spectrogram(signal[:,0], fs = fs, window =name_window, 
								nperseg=nperseg,noverlap = noverlap,
								nfft=None,return_onesided=onesided)

		# Plot Espectrograma + Periodo de Movimientos ---
		cmap=plt.get_cmap('jet')
		num_intervals = np.shape(indexmov[0])[0]

		FIG, AXS = plt.subplots(2,1, gridspec_kw={'height_ratios':[1,4]})

		# Configuration plot
		for i in range(2):
			AXS[i].set_xlim(0,3500)

		AXS[0].set_ylim(0,1)
		AXS[1].set_ylim(10,100)

		AXS[0].set_xticks([])
		AXS[0].set_yticks([])

		AXS[1].set_xlabel('Time [sec]')
		AXS[1].set_ylabel('Frequency [Hz]')
		AXS[1].set_yscale('log')

		# Graphics
		for i in range(num_intervals):
			interval = np.linspace(indexmov[0][i,0]/fs,indexmov[0][i,1]/fs)
			AXS[0].fill_between(interval,1,color='b')

		p = AXS[1].imshow(10*np.log10(Sxx), cmap = 'jet',interpolation = 'bilinear', vmin=-50, vmax=-30,
			origin='lower', extent=[min(t), max(t), min(f), max(f)])

		FIG.colorbar(p, orientation='horizontal')	
   
		plt.show()

		exit()
		# =============================================


#         %% Spectrogram: Comparation between channels from same file.
         
# %         %Create figure
# %         figure('visible','off')
# % 
# %         for n_ch = 1:length(index_channel{NFile})
# %             i=mod(n_ch-1,2)+1;
# %             j=floor((n_ch-1)/2)+1;
# %             
# %             if i ~= 1 
# %                 set(gca,'yticklabel',{[]})
# %             end
# %             
# %             if j ~= 4 
# %                 set(gca,'xticklabel',{[]})
# %             end            
# %         end
# %                
# %        %colorbar
# %        %print(strcat('Results\Spectrogram_',s,NumberFile{NFile}),'-depsc');       
# %        print(strcat('Results/Spectrogram_',s,NumberFile{NFile}),'-depsc');       
        
# %-------------------------------------------------------------------------%
#          %% Temporal evolution of power in particular frequency band: 
#          % Comparation between channels from same file.
 
# %         % Power for each band - OPTION 1
# % %        power_band = zeros(5,length(time),NShowCh);
# % %        
# % %        for band = 1:5
# % %             index_f_band = find(LimBands(band,1) <= f & f < LimBands(band,2));
# % %             power_reduced = power(index_f_band,:,:);
# % %             power_band(band,:,:) = sum (power_reduced,1);
# % %             %power_band = movmean(power_band,4,2);
# % %        end
# %        
# %        disp('Bandas listas')
# %                 
# %         range_time = [IndexUM_Quiet(1,IndexMaxLSignalQ)/fs,IndexUM_Quiet(2,IndexMaxLSignalQ)/fs];
# %         %range_time = [200,340];
# %         %range_time = [0,T(end,1)];
# % 
# %         %Create figure
# %         %figure('visible','off')
# % 
# %         figure()
# %         
# %         for chn = 1:NShowCh
# %             
# %             signal_chn = [signal(end:-1:1,chn); signal(:,chn); signal(end:-1:1,chn)];
# %             [BPFsignal, BPFamp, signalCropped, f_2, t] =...
# %                  functioNStspectrogramFilterHilbert_v1(signal_chn, ...
# %                                                      SFHcfg, indSettling);
# % 
# %   
# %             subplot(NShowCh,1,chn)
# %             hold on
# %             %plot(time,10*log10(power_band(1,:,chn)'),'-b')
# %             %plot(time,10*log10(power_band(2,:,chn)'),'-r')
# %             %plot(time,power_band(2,:,chn)','-r')
# %             plot(T,signal(:,chn),'-g')
# %             plot(t',BPFsignal(:,2),'-r')            
# %             plot(t',BPFamp(:,2).^2,'-b')
# %             %plot(t',movmean(BPFamp(:,1).^2,5*fs),'-b')
# % 
# %             
# %             hold off
# %             %xlim(range_time)
# %                                  
# %             if chn ~= NShowCh 
# %                 set(gca,'xticklabel',{[]})
# %             end
# %             
# %         end
# %        


# %-------------------------------------------------------------------------%
#         %% Spectrogram, Temporal evolution of frequency band, Movement and Artifacts periods
#         % for particular channel.
        
#         %range_time = [IndexUM_Quiet(1,IndexMaxLSignalQ)/fs,IndexUM_Quiet(2,IndexMaxLSignalQ)/fs];
#         %range_time = [200,340];
#         range_time = [0,T(end,1)];
        
#         for chn = 1:NShowCh
            
#             signal_chn = [signal(end:-1:1,chn); signal(:,chn); signal(end:-1:1,chn)];

#             %Power for each band
#             %[BPFsignal, BPFamp, signalCropped, f_2, t] =...
#             %    function_spectrogramFilterHilbert_v1(signal_chn, ...
#             %                                        SFHcfg, indSettling);
                                                
#             %Create figure         
#             %figure('visible','off')
#             figure()
            
#             %Plot Spectrogram
#             subplot('Position',[0.12 0.4 0.8 0.55])
#             max_power = max(max(power(:,:,chn))); 
#             surf(time,f,10*log10(power(:,:,chn)/max_power),'edgecolor','none')
#             view([0 90])
            
#             xlim(range_time)
#             set(gca,'xticklabel',{[]})

#             ylim([0,320])
#             ylabel('Frequecy (Hz)')
#             %set(gca,'Yscale', 'log')
            
#             caxis([-60,0])
#             c=colorbar;
#             c.Location='northoutside';
#             colormap jet
        
# %             %Plot Power Bands 
# %             subplot('Position',[0.12 0.25 0.8 0.1])
# %             hold on
# %             plot(t',BPFamp(:,3).^2,'color',ColorBand{1},'DisplayName',NameBand{3})
# %             hold off
# %         
# %             xlim(range_time)
# %             set(gca,'xticklabel',{[]})        
# %             
# %             aux_ylim = [BPFamp(:,1).^2; BPFamp(:,2).^2];
# %             
# %             ylim([0 30])
# %             ylabel({'Power', 'bands'})
# %             legend('Orientation','Horizontal','Location','northoutside')
        
#             %Plot Movement Periods + Artifacts Index
#             subplot('Position',[0.12 0.1 0.8 0.1])
            
#             hold on
            
#             for i=1:size(IndxUM_Mov,2)
#                  low_bound = linspace(IndxUM_Mov(1,i)/fs,IndxUM_Mov(2,i)/fs);
#                  area(low_bound,up_bound)
#             end

#             plot(T,0.5*index_artifacts,'color','red')

#             hold off
        
#             xlabel('Time (sec)')
#             xlim(range_time) 
        
#             ylabel({'Periods of','movement'})
#             ylim([0,1])
#             set(gca,'yticklabel',{[]})

        
#             %close
            
#        end

