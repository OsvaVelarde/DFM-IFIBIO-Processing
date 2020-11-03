%% BIOMARKERS IN PARKINSON'S DISEASE (PD)

% Author: Osva M Velarde
% Date: 01/07/2018
% Description: 

%------------------------------------------------------------------------%
%% Default configuration
clear global        %Clear all global variables.
clear all           %Clear all local variables.
close all force     %Close all figures (including the wvtool figures).
format long         %Show numbers in long format.
clc                 %Clear screen.
restoredefaultpath

%% Access path to the functions.

pathFunctions1 = '/mnt/BTE2b/DBS/Toolboxs_MATLAB/buzcode-master/externalPackages/FMAToolbox/IO/';
pathFunctions2 = '/mnt/BTE2b/DBS/Toolboxs_MATLAB/buzcode-master/externalPackages/FMAToolbox/Helpers/';
pathFunctions3 = '/mnt/BTE2b/DBS/Toolboxs_MATLAB/Espectrograma/';
pathFunctions4 = '/mnt/BTE2b/DBS/Toolboxs_MATLAB/Normalization/';
pathFunctions5 = '/mnt/BTE2b/DBS/Toolboxs_MATLAB/Spectrum/';
pathFunctions6 = '/mnt/BTE2b/DBS/Toolboxs_MATLAB/Filtering/';
pathFunctions7 = '/mnt/BTE2b/DBS/Toolboxs_MATLAB/Comodulogram/v4/';
pathFunctions8 = '/mnt/BTE2b/DBS/Toolboxs_MATLAB/TimeLockedPlot/';
pathFunctions9 = '/mnt/BTE2b/DBS/Toolboxs_MATLAB/Figure/';

addpath(pathFunctions1, pathFunctions2, ...
        pathFunctions3, pathFunctions4, ...
        pathFunctions5, pathFunctions6, ...
        pathFunctions7, pathFunctions8, ...
        pathFunctions9, ...
        '-begin');

%-------------------------------------------------------------------------%
%% Information about files for the data sets.

LOCAL_PATH_DATA = '/mnt/BTE2b/DBS/Enero-Julio-2019/Biomarkers/Data_Gustavo/data10d2018/';

% Files and Parameters
State = {'Lesionado/20170418-', 'Sham/20170216-'};
NumberStates = length(State);

NumberFile = {'001','002','003','004'};
NNf = length(NumberFile);

NumberChannels = [32,64,64,64];
MAXChn = max(NumberChannels);

%% Electrodes and channels information.

Probe32 = load('/mnt/BTE2b/DBS/2020/GMdata/data/Electrode_32.dat');
Probe64 = load('/mnt/BTE2b/DBS/2020/GMdata/data/Electrode_64.dat');

Probes = {Probe32, Probe64, Probe64, Probe64};

ShowChannels = {1, 20, 62, 62};
%ShowChannels = {[7 23 5 21 3 19 1 17], [6 14 22 30 38 46 54 62], [3 11 19 27 35 43 51 59], [3 11 19 27 35 43 51 59]};
%ShowChannels = {reshape(Probe32',[1,32]), reshape(Probe64',[1,64]), reshape(Probe64',[1,64]), reshape(Probe64',[1,64])};

REFcfg = struct('method','hard','electrode',0);

%% Bands of interest
LimBands = [1,10;51,99]; %0-10-45
NameBands = {'20-30 Hz','230-290 Hz'};
f0 = mean(LimBands,2);
bW = LimBands(:,2)-LimBands(:,1);
bW(2) = max(2*f0(1),bW(2));

%% Movement vs. Quiet information
up_bound = ones(1,100);             

%% Windows parameters
windowName = 'hann';
fs = 1250;           % Hz
sflag = 'symmetric'; %'symmetric' - 'periodic'
windowParam = struct('name',windowName,'sflag',sflag);

%% Checking filter flags.
CHEKING_FILTER_FLAG = 1;
plotFlag = 0;
            
%% Input parameters to compute the comodulogram.
SegmentTimeWidth_CFC = 5 ;%SegmentTimeWidth;
OverlapPercent = 50;

BPFXcfg = struct('Bw',bW(1),... %Parameters for frequency domain.
                 'freqWindowParam',windowParam,...
                 'conv','circular',... %Parameters for time domain.
                 'causal',0,...
                 'zeropadding',0,...
                 'timeWindowParam',windowParam,...
                 'function','function_FDF');
             
BPFYcfg = struct('Bw', bW(2),... %Parameters for frequency domain.
                 'freqWindowParam',windowParam,...
                 'conv','circular',... %Parameters for time domain.
                 'causal',0,...
                 'zeropadding',0,...
                 'timeWindowParam',windowParam,...
                 'function','function_FDF');             

%Comodulogram. For the description of the parameters see function "function_CFCcfg".            
CFCcfg = struct('fXmin',f0(1), 'fXmax',f0(1),...
                'fXres',0.5, 'fXlookAt','PHASE',... 
                'nX',1, 'BPFXcfg',BPFXcfg,...
                'fYmin',f0(2), 'fYmax',f0(2),...
                'fYres',3, 'fYlookAt','PHASEofAMPLITUDE',...
                'nY',1, 'BPFYcfg',BPFYcfg,...
                'saveBPFsignal',0,...                
                'Nbins',18, 'sameNumberOfCycles',1,...
                'CFCmethod','plv', 'verbose',1,...
                'perMethod','sampleShuffling',... %'FFTphaseShuffling',...
                'Nper',0, 'Nrep',0, 'Pvalue',0.05,...
                'corrMultComp','Bonferroni',... %'pixelBased',...
                'fs',fs);

CFCcfg = function_setCFCcfg_v1(CFCcfg);  %% Set configuration of comodulogram

%% Input parameters to compute the Time Locked Index.

%Time Locked Index (TLI).            
TLIcfg = struct('BPFcfg_LF',CFCcfg.fXcfg.BPFcfg,...
                'BPFcfg_HF',CFCcfg.fYcfg.BPFcfg,...
                'abs_HFSignaltimeLockedHFpeaks',0,...
                'abs_HFSignaltimeLockedLFpeaks',0,...
                'LFphase','peaks',...'troughs',...
                'NT',1,...
                'fs',fs,...
                'plot',0);

disp('Loaded parameters')

%% KLMI method
CFCcfg_KLMI = CFCcfg; CFCcfg_KLMI.CFCmethod = 'entropy-based-MI';
CFCcfg_KLMI.fXcfg.Nbins=18;

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%% Processing
for NFile = 1 %1:length(NumberFile) % Loop in files
    for NSts = 1 %1:NumberStates    % Loop in states: sham - lession 
                
        s = State{NSts};         % State.
        NShowCh = length(ShowChannels{NFile}); % Number of channels.  

        %% Load files
        FILENAME_LFP = strcat(LOCAL_PATH_DATA,s,NumberFile{NFile},'.lfp'); 
        FILENAME_MOV = strcat(LOCAL_PATH_DATA,s,NumberFile{NFile},'-IndicesLFP_MovementPeriods.mat');
        FILENAME_QUIET = strcat(LOCAL_PATH_DATA,s,NumberFile{NFile},'-IndicesLFP_QuietPeriods.mat');

        rawSignal = LoadBinary(char(FILENAME_LFP),...
                    'channels', (1:NumberChannels(NFile))',...
                    'nChannels', NumberChannels(NFile)+1,...
                    'start', 0, 'duration', Inf, 'frequency', fs);
        
        load(FILENAME_MOV); %IndxUM_Mov
        load(FILENAME_QUIET); %IndexUM_Quiet
        
        disp('Readed data')
       
        %% Pre-processing: Remove artifacts

        T = (0:+1:size(rawSignal,1)-1).' / fs; %[sec] - Synthesize the time vector. 
        %[index_artifacts, cleanSignal] = ...
        %                remove_artifacts(T,rawSignal, 6250, 0, 0); % Remove artifacts.
        
        %signal = double(cleanSignal);
        signal = double(rawSignal);
        clearvars rawSignal cleanSignal

        %Compute the number of samples of the input signal.
        Nsamples = size(signal,1);
                
        %% Referecing
        REFcfg.electrode = Probes{NFile};
        signal = function_referencing(signal,REFcfg);

        %% Selection of channels for show. 
        signal = signal(:,ShowChannels{NFile});
        
        %% Normalization 
        signal = function_zscore_v1(signal);

        disp('Pre-processing - OK')
        %-------------------------------------------------------------------------%        


        %% Longest signal in Quiet Periods.
        [MaxLenSignalQ, IndexMaxLSignalQ] = max(IndexUM_Quiet(2,:)-IndexUM_Quiet(1,:));
 
        % ----------------------------------------------------------------------- %                                       
        %% Check the filters and compute the settling time (percLevel).

        %Initialize the index.
        indSettling = NaN;
            
        if CHEKING_FILTER_FLAG == 1, %Check the filters just once.

            CHEKING_FILTER_FLAG = 0;
            
            %BPFs for the "x" axis of the comodulogram.
            indSettlingLF = function_checkFilter_v1(CFCcfg.fXcfg.BPFcfg, CFCcfg.fs, Nsamples, plotFlag);           
            
            %BPFs for the "y" axis of the comodulogram.
            indSettlingHF = function_checkFilter_v1(CFCcfg.fYcfg.BPFcfg, CFCcfg.fs, Nsamples, plotFlag);
            
            %Compute the maximum settling time.
            indSettling = max(max([indSettlingHF; indSettlingLF.']));
            
        end %Check the filters just once.

        %Compute the maximum settling time.
        indSettling = max([indSettling, Nsamples+1]);

        disp('Check Filters - OK')

%-------------------------------------------------------------------------%
       %% Loops across the channels
        for Nch = 1:NShowCh  % Loop across the channels.
           disp(['Canal ',int2str(Nch)])

           SignalCh = signal(:,Nch);            

           %% Compute the Band-Pass Filtering.

           % Reflect the time series to minimize edge artifacts due to the transient response of the BPFs.
           SignalCh = [SignalCh(end:-1:1); SignalCh; SignalCh(end:-1:1)];
                  
           [LFSIGNAL, ~] = function_comodulogramBPF_v1(SignalCh, ...
                               CFCcfg.fXcfg.BPFcfg, CFCcfg.fs, indSettling);
                           
           [HFSIGNAL, ~] = function_comodulogramBPF_v1(SignalCh, ...
                               CFCcfg.fYcfg.BPFcfg, CFCcfg.fs, indSettling);
                      
           % Restore the length of the raw signal.
           SignalCh = SignalCh(indSettling:end-(indSettling-1),:);
           
          
           disp('LF and HF signal - OK')
           
           %% Segmentation of the Band-Pass Filtered time series.
           indSegment = function_segmentation_v1(SegmentTimeWidth_CFC, OverlapPercent, Nsamples, fs);

           NumSegments = size(indSegment,1);
           
           IndexesEvolution = zeros(NumSegments,6);           
                             
%           for ii=1:+1:size(indSegment,1), %Loop across the segments.
           for ii=1:2, %Loop across the segments.
                
                disp(['Processing segment ',num2str(ii),' out of ',num2str(NumSegments),'.'])
               
               %Non-causal time vector
               taux2 = mean(T(indSegment(ii,1):indSegment(ii,2)));

               %Compute the Raw segment -----------------------------------
               segRawSignal = SignalCh(indSegment(ii,1):indSegment(ii,2));

               %Compute the LF segment ------------------------------------
               segLFSIGNAL = LFSIGNAL(indSegment(ii,1):indSegment(ii,2),:);

               %Compute the HF segment ------------------------------------
               segHFSIGNAL = HFSIGNAL(indSegment(ii,1):indSegment(ii,2),:,:);
 
%                %Compute the TLI -------------------------------------------
% 
%                TLI = NaN(length(CFCcfg.fYcfg.BPFcfg.f0),length(CFCcfg.fYcfg.BPFcfg.Bw));
%                for indf=1:+1:length(CFCcfg.fYcfg.BPFcfg.f0), %Loop for Frequencies.
%                    for indBw=1:+1:length(CFCcfg.fYcfg.BPFcfg.Bw), %Loop for Bandwidths.
%                        
%                        %TLI
%                        TLIcfg_local = TLIcfg;
%                        TLIcfg_local.BPFcfg_HF.f0 = TLIcfg.BPFcfg_HF.f0(indf);
%                        TLIcfg_local.BPFcfg_LF.f0 = TLIcfg.BPFcfg_LF.f0(indBw);
%                        TLIcfg_local.BPFcfg_HF.Bw = TLIcfg.BPFcfg_HF.Bw(indBw);
%     
%                        [rawSignaltimeLokedHFpeaks, rawSignaltimeLokedLFpeaks,...
%                         HFSignaltimeLokedHFpeaks, HFSignaltimeLokedLFpeaks,...
%                         LFSignaltimeLokedHFpeaks, LFSignaltimeLokedLFpeaks,...
%                         TLI(indf,indBw), indPeak_HF, indPeak_LF, sampleT] =...
%                         function_TimeLockedIndex_v1(segRawSignal, segLFSIGNAL(:,indBw), squeeze(segHFSIGNAL(:,indf,indBw)), TLIcfg_local);
%     
%                    end
%                end
%                TLI = squeeze(TLI);

               %-------------------------------------------------------------------------%
               %% Compute features ---------------------------------------

               %Compute the phase/amplitude/frequency for the LF segment -----------------
               %Compute the number of samples of the segment.
               Namples_per_segment = size(segLFSIGNAL,1);

               %Z-score normalization of the segment (in order to have zero mean and unit variance).
               segLFSIGNAL = function_zscore_v1(segLFSIGNAL);

               %Reflect the segment to minimize edge artifacts due to the transient response of the Hilbert transform.
               segLFSIGNAL = [segLFSIGNAL(end:-1:1,:); segLFSIGNAL; segLFSIGNAL(end:-1:1,:)];

               %Compute the phase/amplitude/frequency segment.
               [segX, ~] = function_comodulogramFeature_v1(segLFSIGNAL, CFCcfg.fXcfg, CFCcfg.fs, Namples_per_segment+1);

               %Compute the amplitudes envelope for the weighted PLV.
               fXcfg_aux = CFCcfg.fXcfg; fXcfg_aux.lookAt = 'AMPLITUDE';
               LFamp = function_comodulogramFeature_v1(segLFSIGNAL, fXcfg_aux, CFCcfg.fs, Namples_per_segment+1);
               clear fXcfg_aux

               %--------------------------------------------------------------------------

               %Compute the phase/amplitude/frequency for the HF segment -----------------
               %Z-score normalization of the segment (in order to have zero mean and unit variance).
               for jj=1:+1:length(CFCcfg.fYcfg.BPFcfg.Bw); %Loop for Bandwidths.
                    segHFSIGNAL(:,:,jj) = function_zscore_v1(segHFSIGNAL(:,:,jj));
               end

               %Reflect the segment to minimize edge artifacts due to the transient response of the Hilbert transform.
               segHFSIGNAL = [segHFSIGNAL(end:-1:1,:,:); segHFSIGNAL; segHFSIGNAL(end:-1:1,:,:)];

               %Compute the phase/amplitude/frequency segment.
               [segY, ~] = function_comodulogramFeature_v1(segHFSIGNAL, CFCcfg.fYcfg, CFCcfg.fs, Namples_per_segment+1);

               %Compute the amplitudes envelope for the weighted PLV.
               fYcfg_aux = CFCcfg.fYcfg; fYcfg_aux.lookAt = 'AMPLITUDE';
               HFamp = function_comodulogramFeature_v1(segHFSIGNAL, fYcfg_aux, CFCcfg.fs, Namples_per_segment+1);
               clear fYcfg_aux
               
               %--------------------------------------------------------------------------
               %% Save Total segments
               %TotalSegX = [TotalSegX;segX];
               %TotalSegY = [TotalSegY;segY];

               %figure()
               %plot(segHFSIGNAL)
               
               %figure()
               %plot(segY)
               %% Compute the Phase Locking Value.
               [PLV, wLF_PLV, wHF_PLV] = function_PLV_v1(segX, segY, [], [], CFCcfg);
               
%                if CFCcfg.Nrep && CFCcfg.Nrep, 
% 
%                    % Compute surrogate Phase Locking Values via single-trial non-parametric permutation methods.
%                    %Memory pre-allocation for speed up the loop.
%                    PLVs = NaN(length(CFCcfg.fYcfg.BPFcfg.f0),length(CFCcfg.fXcfg.BPFcfg.f0),CFCcfg.Nper*CFCcfg.Nrep);
%                    wLF_PLVs = PLVs;
%                    wHF_PLVs = PLVs;
%                    for jj=1:+1:CFCcfg.Nrep, %Loop across the repetitions.
%                         [PLVs(:,:,(jj-1)*CFCcfg.Nper+1:jj*CFCcfg.Nper),...
%                          wLF_PLVs(:,:,(jj-1)*CFCcfg.Nper+1:jj*CFCcfg.Nper),...
%                          wHF_PLVs(:,:,(jj-1)*CFCcfg.Nper+1:jj*CFCcfg.Nper)] =...
%                         function_comodulogramPermutation_v1(segX, segY, LFamp, HFamp, CFCcfg);
%                    end
% 
%                    % 9) Surrogate Control Analysis.
%                    [absPLVz, absPLVth, absPLVs_mean, absPLVs_std, absPLVs_threshold] = function_comodulogramSCA_v1(abs(PLV), abs(PLVs), CFCcfg);
%                    [abs_wLF_PLVz, abs_wLF_PLVth, abs_wLF_PLVs_mean, abs_wLF_PLVs_std, abs_wLF_PLVs_threshold] = function_comodulogramSCA_v1(abs(wLF_PLV), abs(wLF_PLVs), CFCcfg);
%                    [abs_wHF_PLVz, abs_wHF_PLVth, abs_wHF_PLVs_mean, abs_wHF_PLVs_std, abs_wHF_PLVs_threshold] = function_comodulogramSCA_v1(abs(wHF_PLV), abs(wHF_PLVs), CFCcfg);
% 
%                else
% 
%                        PLVs = NaN;
%                        absPLVz = NaN;
%                        absPLVth = NaN;
%                        absPLVs_mean = NaN;
%                        absPLVs_std = NaN;
%                        absPLVs_threshold = NaN;
%                        wLF_PLVs = NaN;
%                        abs_wLF_PLVz = NaN;
%                        abs_wLF_PLVth = NaN;
%                        abs_wLF_PLVs_mean = NaN;
%                        abs_wLF_PLVs_std = NaN;
%                        abs_wLF_PLVs_threshold = NaN;
%                        wHF_PLVs = NaN;
%                        abs_wHF_PLVz = NaN;
%                        abs_wHF_PLVth = NaN;
%                        abs_wHF_PLVs_mean = NaN;
%                        abs_wHF_PLVs_std = NaN;
%                        abs_wHF_PLVs_threshold = NaN;
% 
%                end %if CFCcfg.Nrep && CFCcfg.Nrep,
               
               %% KLMI
               %p = function_comodulogramHistogram_v1(segX, HFamp, CFCcfg_KLMI);
               %KLMI = function_KLMI_v1(p, CFCcfg_KLMI);

               IndexesEvolution(ii,1)=taux2;
               IndexesEvolution(ii,2)=abs(PLV);
               %IndexesEvolution(ii,2)=TLI;
               %IndexesEvolution(ii,3)=KLMI;
               
            end %Loop across the segments.

            %-----------------------Plots and Results --------------------------------%      

            %% Spectrogram + HF amplitude + PLV + TLI + Movement and Artifacts periods
            % for particular channel.
        
            %rangeTime = [IndexUM_Quiet(1,IndexMaxLSignalQ)/fs,IndexUM_Quiet(2,IndexMaxLSignalQ)/fs];
            %rangeTime = [515,530];
            rangeTime = [0,T(end,1)];
                                                            
            % Create figure         
            %figure('visible','off')
            figure()
            
            % Plot Spectrogram --------------------------------------------
            subplot('Position',[0.12 0.55 0.8 0.45])
            %max_power = max(max(power(:,:,Nch))); 
            %surf(time,f,10*log10(power(:,:,Nch)/max_power),'edgecolor','none')
            %view([0 90])
            
            xlim(rangeTime)
            set(gca,'xticklabel',{[]})
        
            ylim([0,300])
            ylabel('Frequecy (Hz)')
            set(gca,'Yscale', 'log')
            
            caxis([-60,0])
            c=colorbar;
            c.Location='northoutside';
            colormap jet
        
            % Plot PLV - TLI ----------------------------------------------
            subplot('Position',[0.12 0.25 0.8 0.1])
            hold on
            plot(IndexesEvolution(:,1),IndexesEvolution(:,2),'DisplayName','PLV')
            %plot(IndexesEvolution(:,1),IndexesEvolution(:,5),'DisplayName','TLI')
            hold off
            
            xlim(rangeTime)
            set(gca,'xticklabel',{[]})        
                        
            %ylim([0 0.4])
            ylabel({'Indexes'})
            legend('Orientation','Horizontal','Location','north')
            
            % Plot Movement Periods + Artifacts Index ---------------------
            subplot('Position',[0.12 0.1 0.8 0.1])
            hold on
            
            for i=1:size(IndxUM_Mov,2)
                 low_bound = linspace(IndxUM_Mov(1,i)/fs,IndxUM_Mov(2,i)/fs);
                 area(low_bound,up_bound)
            end

            plot(T,0.5*index_artifacts,'color','red')

            hold off
        
            xlabel('Time (sec)')
            xlim(rangeTime) 
        
            ylabel({'Periods of','movement'})
            ylim([0,1])
            set(gca,'yticklabel',{[]})

            % Save figure (Windows and Ubuntu VERSIONs)
            %savefig(strcat('Results/Spectrogram_',s,NumberFile{NFile},'_Channel_',int2str(ShowChannels{NFile}(chn)),'.fig'))

            %print('Results/Prueba','-depsc');
        
            %close
            %-------------------------------------------------------------------------%

       end 
    end
end