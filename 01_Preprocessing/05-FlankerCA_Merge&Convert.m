function FlankerCA_Merge&Convert
% NOTES:
% "whatever exp name"_TimeFreqAnalysis is used for EEG data with bin of
% "ALL", in which power and phase of all epochs are extracted altogether.
% Analyze data with conditionining based on the behavioral data or factors
% that were not coded in the event codes!!!


% TO DO:
% 1.Merge EEG and Behavioral data sets for each subject
% 2.Convert to H5 format
% *************************************************************************

global EXCON FEMODEL CLASSIFY MI FILEPATH ENTROPY

%%
clc
clear all
close all

% pool = parpool(2)

% run all prep (includes eeglab stuff and time-frequency decomposition
%PrepAll_EEGLAB_BP_Analysis_MM % NOTE: masterdir should be 'EEG'

% Basic Setting
EXCON.ANALYSIS='WAVELET';
EXCON.NORMALIZE='dB';
EXCON.EXCODE='FlankerConj';
EXCON.EPOCHDUR={-200 1200};%feedbacl(700)+trial(2000)->first 1000 ms after stim onset is included!
EXCON.BLOCKCODE=[200,-7];%these codes should be dependent on the ev codes!
EXCON.TRIALCODE=[100,-8];%it is identifying trial from ev code for EEG data.
EXCON.SUBSKILLED=[];%
EXCON.STIMLOCKED=false;

% Forward Endocing Model Setting for task and serial position
% FEMODEL.SRATE=250;
% FEMODEL.CLASSWINDOW=1;
% FEMODEL.PHASEWINDOW=10;
% FEMODEL.PHASETIMER=[];
% CLASSIFY.SRATE=250;
% CLASSIFY.CLASSWINDOW=25;
% CLASSIFY.PHASEWINDOW=2;
% CLASSIFY.PHASETIMER=[350,450];
% MI.SRATE=250;
% MI.CLASSWINDOW=100;

% Path Setting
FILEPATH.MASTERDIR=uigetdir; %NOTE: should be 'DataAnalysis/Matlab'
FILEPATH.BEHDIR=strcat(FILEPATH.MASTERDIR,'/BEH/w_ALLGRAND');
FILEPATH.EEGDIR=strcat(FILEPATH.MASTERDIR,'/EEG');
FILEPATH.ALLGRANDDIR=strcat(FILEPATH.EEGDIR,'/w_ALLGRAND');
cd(FILEPATH.EEGDIR);AllSubInfo=dir('A*');
for i=1:length(AllSubInfo),
    subs(i)=str2num(strrep(AllSubInfo(i).name,'A',''));
    FILEPATH.SUBDIR(i)={strcat(FILEPATH.EEGDIR,'/',AllSubInfo(i).name)};
end

% Load up  BEH data
cd(FILEPATH.BEHDIR)
ALLDATA=v2struct(load('processedDataS.mat'));
ALLDATA.raweeg(:,1)=cell(size(ALLDATA,1),1);
ALLDATA.eegpower(:,1)=cell(size(ALLDATA,1),1);
ALLDATA.eegphase(:,1)=cell(size(ALLDATA,1),1);
ALLDATA.reject(:,1)=nan(size(ALLDATA,1),1);
ALLDATA.SUBID=ALLDATA.ID;

ALLDATA=ALLDATA(ismember(ALLDATA.Practice,0),:);%double check to remove pracitice blocks 
% % Entropy
% ENTROPY.toi=[0.375];%all across the trials
% ENTROPY.hwinsize = 0.375;
% ENTROPY.m=2;% pattern length  commonly set to 2
% ENTROPY.r= 0.5; % similarity criterion ommonly ranges between .1 and .5
% ENTROPY.scales= [20:42]; % scale list (e.g.,1:42)
% ENTROPY.depV={'eegraw'};
% cd(FILEPATH.EEGDIR)
% ENTROPY.EEGt=v2struct(load('FTS_fieldtrip_temp.mat'));%template

% Coding group
subsKeep=~ismember(subs,EXCON.SUBSKILLED);
FILEPATH.SUBDIR=FILEPATH.SUBDIR(subs);%FILEPATH.SUBDIR=FILEPATH.SUBDIR(subsKeep);
subs=subs(subsKeep);
EXCON.SUBS=subs;


%%
% =========================================================================
% Prepare template fieldtrip EEG (for mse_core)
% =========================================================================

% % % To run mse_core, make a file formatted in fildtrip style
% % % Load this template, and adjust trial and time field!
% cd('whoever subject')
% EEGf=v2struct(load('A200_EEG_preprocessed.mat'));
% EEGf=eeglab2fieldtrip(EEGf,'preprocessing');
% EEGf.trial={nan(size(EEGf.trial{1}))};
% EEGf.time=EEGf.time(1);
% cd(FILEPATH.EEGDIR)
% saveMatFile(EEGf,'EEGf','FTS_fieldtrip_temp');



%%
% =========================================================================
% PROCESS 2 and 3: Calculate all sorts of things! Currently, functions are
% available for mean power, inter-trial phase clustering,inter-time phase
% clustering, bivariate phaseband coherence, etc....
% Also tf_sorcontraipsi is available for sorting lateralized EEG data!
% =========================================================================
clc
tic

% open parallel port if it is not open
paralalleStatus=gcp;
if ~paralalleStatus.Connected,parpool(3);end

% recency bias is z 1 ! all others are z2
%ct=1;

for sub=1:length(subs)
    
    s=subs(sub);
    sprintf('starting subject %d',s) 

    % =====================================================================
    % PROCESS 2: Make sure template data exist in workspace...
    % =====================================================================
    % Loading up data and template ----------------------------------------
    %Specify directory   
    FILEPATH.SUBJECT=uniformLabel({'A','0',3},s,'first');
    FILEPATH.DATASETDIR=strcat(FILEPATH.EEGDIR,'/',FILEPATH.SUBJECT{:},'/DATASETS');
    %FILEPATH.SUBJECT=FILEPATH.SUBDIR{ct};
    %FILEPATH.DATASETDIR=strcat(FILEPATH.SUBJECT,'/DATASETS');
    data2load={'eegpower'};%{'eegraw','eegpower','eegphase'}
    
%     if s > 6
%         % Block and trial codes are switched after subject 7
%         EXCON.BLOCKCODE=[100,-6];%these codes should be dependent on the ev codes!
%         EXCON.TRIALCODE=[200,-5];%it is identifying trial from ev code for EEG data.
%     end
    
    [DATA,EEG]=tf_mergeEEGBEH_Conj_MM(s,ALLDATA,FILEPATH,EXCON,data2load);%for ERP
    
    % =====================================================================
    % Filter data
    % =====================================================================
    DATA=DATA(DATA.Block>1,:);
    DATA=DATA(DATA.Error ~= 1 & DATA.L1Error ~= 1 & DATA.RT<=2000,:);
    DATA=DATA(DATA.Trial>1,:);
    % =====================================================================
    % Regress data????
    % =====================================================================
    

    % =====================================================================`
    % 
    % Basline data???? 
    % =====================================================================
    
    %Change setting
    %switch data2load{:}
    %  case {'eegraw'},ntype='subtract';bInd=dsearchn(EEG.BINEPOCH.times',[-200,0]')';         
    %  case {'eegpower'},ntype='dB';bInd=dsearchn(EEG.BINEPOCH.times',[-350,-150]')';    
    %end
    
    %Normalization
    %DATA.(data2load{:})=tf_basenormalize(DATA.(data2load{:}),ntype,bInd);

    
    % =====================================================================
    % PROCESS 2: Prep data to transfer to R!!
    % =====================================================================
    %tf_preptrans2R(DATA,EEG,FILEPATH,'depV',data2load,'timeC',[-200,1500],'electrpde',{'ALL'});
    

    %Setting
    SET.subInd=find(subs==s);
    SET.depV={'eegpower'};
    SET.fband={'thetaFreqs','alphaFreqs','betaFreqs','gammaFreqs'};
    SET.elec={'ALL'};
    SET.zscore=true;
    SET.fext={'SINGLE_allFB'};
    tf_trans2env_MM(DATA,EEG,FILEPATH,SET);

end
