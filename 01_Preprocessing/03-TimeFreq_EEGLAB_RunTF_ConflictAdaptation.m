function TimeFreq_EEGLAB_RunTF_ConflictAdaptation(PREP,varargin)
% NOTES:
% 1, freq band is updated!!! reanalyze the following experiments!!
% (switchload,updatefrontal,repeatencode,carryover)
% PROCESS TOWARD OSCILLATORY ANALYSIS!!------------------------------------
% % -----------------------------------------------------------------------
% STEP1:
% Complete preprocessing, binning, and artifact rejection through EEGLAB
% Make sure to save EEG file without any filtering (baselining should not influence)
% 
% STEP2:
% Create summary of epochs specific to a unique bin. Copy all of useful and
% necessary information from EEG feilds from STEP1
% 
% STEP3:
% Extract amplitude(power), phase of EEG signal based on either 1)complex
% Morlet wavelet analysis, 2)Bandpath-Hilbert transform, 3)short-time Fast 
% Fourier transform, for all channels and all trials. Get average.

% *************************************************************************

%%
fprintf('~~~~~Starting TF analysis for %s !!~~~~~\n.',PREP.SUBLABEL{PREP.i});

% Make sure to have EEG data either from input or from mat file!
cd(PREP.SUBDIR{PREP.i});if ~exist(fullfile(cd, 'DATASETS'),'dir'),mkdir('DATASETS');end;PREP.DATASETDIR=strcat(PREP.SUBDIR{PREP.i},'/DATASETS');
cd(PREP.DATASETDIR);if isempty(varargin),EEG=v2struct(load([PREP.SUBLABEL{PREP.i},'_EEG_preprocessed.mat']));end;

%%
% Transfer basic information about data!!To be consistent with the other functions,
% the file structure is based off of the older version of EEGLAB...!!
% Some stuff are redundant, but better than adjusting everything!

% Setting 
% EEG.BINEPOCH.binname=actBin.description;
% EEG.BINEPOCH.bincode=strcat(actBin.namebin(1),actBin.namebin(end));
EEG.BINEPOCH.numepochs= EEG.trials;
EEG.BINEPOCH.pnts= EEG.pnts;
EEG.BINEPOCH.srate=EEG.srate;
EEG.BINEPOCH.chanlocs= EEG.chanlocs;
EEG.BINEPOCH.rateAcq=1000/EEG.srate;
EEG.BINEPOCH.times=EEG.times;
EEG.BINEPOCH.preDP=length(find(EEG.times<0));
EEG.BINEPOCH.postDP=length(find(EEG.times>=0));
EEG.BINEPOCH.allFreqs=PREP.FREQINFO.allFreqs;
EEG.BINEPOCH.numFreqs=PREP.FREQINFO.numFreqs;
% EEG.BINEPOCH.normalize=EEG.EXCON.NORMALIZE;
% EEG.BINEPOCH.baseline=EEG.EXCON.BASELINE;

% EEG data and bins
bind=true(size(EEG.epoch,2),1);%grab everything in a bin!
EEG.BINEPOCH.trials=length(find(bind));
EEG.BINEPOCH.bepoch=EEG.epoch(bind);
EEG.BINEPOCH.eegraw=squeeze(EEG.data(:,:,bind));%too huge!
EEG.BINEPOCH.reject=EEG.reject.rejmanual(1,bind);
% bind= false(size(EEG.epoch,2),1);for e=1:size(EEG.epoch,2),bind(e,:)=any(find(cell2mat(EEG.epoch(:,e).eventbini)==b));end

%keyboard;
%%
%TIME-FREQUENCY ANALYSIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
switch PREP.TFANALYSIS
    case {'WAVELET'},EEG.BINEPOCH.(PREP.TFANALYSIS)=tf_cmwavelet_n2(EEG.BINEPOCH,'wvgsd',[3,10]);%tf_cmwavelet;old function!
    %case {'HILBERT'},%simply apply hilbert after filtering of concatanated data
    %case {'sFFT'},%do FFT for each epoch with sliding window
end

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%Pros and Cons of WAV?ELET, HILBERT, and sFFT
%1,WAVLET
%pros=fast, smoothed, control over time & frequency precision,suited
%for non-stationary transient oscillation
%cons=no control over filtering kernel
%2,BANDPATH and HILBERT
%pros= could use a variety of filters to band-path signal
%cons=Matlab hilbert and filtering is slow
%3,SHORT TIME FFT
%pros= no edge artifact
%cons=slowest among all
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%%
%To make a file smaller!!Keep only absolutely necessary....!!
EEG_r=struct;
% EEG_r.EXCON=EEG.EXCON;
% EEG_r.FILEPATH=EEG.FILEPATH;
% EEG_r.BININFO=EEG.BININFO;
EEG_r.EVENTLIST=EEG.EVENTLIST;
EEG_r.FREQINFO=PREP.FREQINFO;
EEG_r.BINEPOCH=EEG.BINEPOCH;
EEG_r.BINEPOCH.bepoch=[];

% % Save EEG , POWER and PHASE into bin/h5 file!
fname=strcat(PREP.DATASETDIR,filesep,strcat(PREP.SUBLABEL{PREP.i},'_EEG_set'));
tf_savebins(fname,EEG_r,{'eegraw','eegpower'},'tepochdur',[PREP.TEPOCHDUR{:}]);
% tf_savebins(fname,EEG_r,{'eegraw'},'tepochdur',[PREP.TEPOCHDUR{:}]);
% tf_saveh5(fname,EEG_r,{'eegraw'},'tepochdur',[PREP.TEPOCHDUR{:}]);

fprintf('~~~~~Finished TF analysis for %s !!~~~~~\n.',PREP.SUBLABEL{PREP.i});

end