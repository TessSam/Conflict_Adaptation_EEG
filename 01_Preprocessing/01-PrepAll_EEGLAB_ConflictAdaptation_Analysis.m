function PrepAll_EEGLAB_ConflictAdaptation_Analysis

%%
% ccc
clc
clear all
close all


% Settings
% PREPROCESS SETTING-------------------------------------------------------
% -------------------------------------------------------------------------
% Experiment & Path 
addpath(genpath(pwd)); % make sure this command happens while you're in the Dropbox (need to add lib functions and all functions from P.FlankerConj)
PREP=struct;

PREP.MASTERDIR=uigetdir; %set to DataAnalysis>Matlab>EEG
%PREP.MASTERDIR='/Volumes/Data2/ulrichmayr/Dropbox/P.FlankerConj';%/DataAnalysis/Matlab/EEG';% for Tess CP
prepDir=which('TimeFreq_EEGLAB_Preprocess_Conj_MM.m');
cd(PREP.MASTERDIR);allsubs=dir('A*');
addpath(genpath(pwd));
%PREP.SUBDIR=[PREP.MASTERDIR '/eeg_data/A005'];
PREP.SUBLABEL={allsubs.name};
PREP.SUBID=cell2mat(cellfun(@(sd)str2double(sd(regexp(sd,'\d'))),{allsubs.name},'Uni',false));
PREP.SUBDIR=cellfun(@(sd)strcat(PREP.MASTERDIR,filesep,sd),{allsubs.name},'Uni',false);

% Preprocessing steps
% PREP.DOCALIB=true;%calibration (oregon specific)
PREP.DOREREF=false;%A1-A2 mastoid reference adjustment
PREP.DORESAMPLE=true;% re(down)sampling
PREP.STIMLOCKED=false;% if true, epoching is stimulus-locked (flanker onset), otherwise, response-locked
PREP.DOEPOCH=true;%epocing specifedd by epoch duration
PREP.DOBIN=true;%binning based on BDF
PREP.DOBLOCK=true;%check blocking (no AR!)
PREP.DOEYEMOVE=true;%eye movement detection(no AR!)
PREP.DOEYEBLINK=true;%blink detection(no AR!)
PREP.DOFILTER=true;%filtering of signal
PREP.SAVESET=false;%save result into .set (EEGLAB files) 

% Epocing and Binning
if PREP.STIMLOCKED
    PREP.EPOCHDUR={-200 1200};%to be epoched range of data % {-200 1200} for flanker stimulus-locked | {-1000 400} for response-locked
    PREP.AFEPOCHDUR={0 1200};%to be checked for artifact rejection % {0 1200} for stimulus-locked | {-800 200} for response-locked
    PREP.TEPOCHDUR={-200 1200};%to be stored range of data after TF analysis % {-200 1200} for stimulus-locked | {-1000 400} for response-locked
    PREP.BDFNAME='BDF_A.txt';%list all BDFs to use % 'BDF_A.txt' for stimulus-locked | 'BDF_R.txt' for response-locked
else
    PREP.EPOCHDUR={-800 200};%to be epoched range of data % {-200 1200} for flanker stimulus-locked | {-1000 400} for response-locked
    PREP.AFEPOCHDUR={-800 0};%to be checked for artifact rejection % {0 1200} for stimulus-locked | {-800 200} for response-locked
    PREP.TEPOCHDUR={-800 200};%to be stored range of data after TF analysis % {-200 1200} for stimulus-locked | {-1000 400} for response-locked
    PREP.BDFNAME='BDF_R.txt';%list all BDFs to use % 'BDF_A.txt' for stimulus-locked | 'BDF_R.txt' for response-locked
    PREP.SUBLABEL=strcat(PREP.SUBLABEL,'_R');
end
PREP.SRATE=1000;
PREP.BASELINE='none';%none,pre
PREP.BDFDIR=strcat(PREP.MASTERDIR,filesep,PREP.BDFNAME);%location of BDF.txt file
% PREP.BOPNAME='BOP.txt';%list all BDFs to use
% PREP.BOPDIR=strcat(PREP.MASTERDIR,filesep,PREP.BOPNAME);%location of BOP.txt file
% BDF = bin description file summarising bin(s), which is group of epochs 
% BOP = bin operation file specificing what to do with bins (not used!)
% Both files should be adjusted for specific experiment!!!!!!!!!!!!!!!!!!!

% BP 32ch system 
%PREP.ELECDIR=which('AS-32_NO_REF.bvef');%location of channel info file
%PREP.ELECDIR=[PREP.MASTERDIR filesep 'AS-32_NO_REF.bvef'];
%PREP.ELECDIR='AS-32_NO_REF.bvef';
PREP.ELECDIR='CMA-32_REF.bvef';
PREP.ELECREMAP={{'TP9','A1'},{'TP10','A2'}}; % Name of workstation!
PREP.ELECREFS={'A1','A2'};

% Re(down)sampling
PREP.RESAMPLE=250;% in Hz (originally 250; 1 sample every 4ms)

% Blocking-Flat line
PREP.FLATCHAN=[2:30];%electrodes to act on
PREP.FLATWIND=200;%window size (ms)
PREP.FLATCUT=[-0.1,0.1];%threshold of EEG value to be excluded (mV)
PREP.FLATFLAG=[1,2];%flag code used in artifact rejection summary

% Blink 
PREP.BLINKCHAN=[1];%electrodes to act on
PREP.BLINKWIND=200;%window size (ms)
PREP.BLINKSTEP=50;%stepping window size (ms)
PREP.BLINKCUT=[120];%threshold of EEG value to be excluded (mV) (80~200); originally set to 125
PREP.BLINKFLAG=[1,3];%flag code used in artifact rejection summary

% Eye movemenet 
PREP.EYEMCHAN=[31];%electrodes to act on
PREP.EYEMWIND=200;%window size (ms)
PREP.EYEMSTEP=10;%stepping window size (ms)
PREP.EYEMCUT=[75];%threshold of EEG value to be excluded (mV) (50~100); originally set to 75
PREP.EYEMFLAG=[1,4];%flag code used in artifact rejection summary

% Filtering
PREP.FILTERCHAN=[2:30];%electrodes to act on
PREP.FILTERCUT=[0.01,45];%to be included range of signals 



% TIME FREQUENCY ANALYSIS SETTING------------------------------------------
% -------------------------------------------------------------------------
FREQ.deltaFreqs=[1:3];
FREQ.thetaFreqs=[4:7];
FREQ.alphaFreqs=[7:12];
FREQ.betaFreqs=[13:31];
FREQ.gammaFreqs=[31:40];
FREQ.allFreqs=[2:1:40];
FREQ.numFreqs=[length(FREQ.allFreqs)];
PREP.FREQINFO=FREQ;
PREP.TFANALYSIS='WAVELET';
% Depending on exact method, default settings are defined later....!
% It would be nice to adjsut the default from here!!

%%
% LOOP FOR EACH SUBJECT----------------------------------------------------
% -------------------------------------------------------------------------
%clc

for sub=1:length(PREP.SUBID)
%for sub=1

    tic
    
    cd(PREP.MASTERDIR)
    %Update loop specific setting
    PREP.s=PREP.SUBID(sub);PREP.i=sub;
    
    %1,Run preprocessing and extract EEG
    [EEG]=TimeFreq_EEGLAB_Preprocess_ConflictAdaptation(PREP);
    
    %2,Time frequency analysis
    TimeFreq_EEGLAB_RunTF_ConflictAdaptation(PREP);
    
    toc
end

close all
disp('~~~~~All done!!!!~~~~~')











end