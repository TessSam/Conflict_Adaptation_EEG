function FlankerCA_BehAnalysis
%==========================================================================
% UPDATES: Before running this script, use terminal to complete the following steps:
% 1. cd into beh data folder ('/Users/Melissa/Dropbox (University of Oregon)/P.FlankerConj/DataAnalysis/Matlab/BEH')
%
% 2. concatenate all behavioral data (text files) using the following
% command: awk 'FNR != 1' *.txt > FlankerConj_Beh.txt

% 3. add the variable names back in as first row

% 4. copy the file into w_ALLGRAND as FlankerConj_Beh.txt
%==========================================================================
% Or
% cd into beh data folder and run the following terminal commands
% awk 'FNR != 1' *.txt > ConflictFlanker_allBehData.txt  
% awk '!($9="")' FlankerConj_allBehData.txt > FlankerConj_Beh.txt

% NOTES:
% 
%==========================================================================
% BASIC STEPS--------------------------------------------------------------
% 1, Join response related data and parameters into one dataset
% ->behavioralDS

% 2, Prepare additional parameters
% ->refinedDS,refinedDS_1,refinedDS_2 and etc

% 3, Apply cut off 
% ->processedDS

% 4, Do statistical analysis and Plotting
% (if needed, create a data spreadsheet for R analysis!!)


%%

%BEHAVIORAL DATA-----------------------------------------------------------
%Analysis Start!===========================================================
%==========================================================================



%%
clc
clear all
close all

addpath(genpath('~/Dropbox (University of Oregon)/libMatlabFunctions_andR'));

% Coding Basic Info
BASIC.MASTERDIR=uigetdir; %go into behavioral data folder (in DataAnalysis/BEH here)
BASIC.EXCODE='cf';
BASIC.FILETOLOAD_R='FlankerConj_Beh';
BASIC.ALLGRANDIR=strcat(BASIC.MASTERDIR,'/w_ALLGRAND');

% Get all subject directories...
cd(BASIC.ALLGRANDIR);
%BASIC.PARAMFILES=dir('*PARAM_ST.mat');%first 3 digits are subjec ID
%subs=cellfun(@(f)str2double(f(idx(regexp(f,'\d'),[1:3]))),{BASIC.PARAMFILES.name},'Uni',false);
%BASIC.ALLSUBS=unique(vertcat(subs{:}));

%%
% Create Beh mat file
ds=readtable(strcat(BASIC.FILETOLOAD_R,'.txt'), 'Delimiter', '\t', 'ReadVariableNames', true);
saveMatFile(ds,'ds','behDataS');
%behDataS=importdata(strcat(BASIC.FILETOLOAD_R,'.txt'), ' ', 1);
%txt2dataSet(strcat(BASIC.FILETOLOAD_R,'.txt'),'behDataS');

%%

% % %Process2:Preparing all variables in dataset array
cd(BASIC.ALLGRANDIR)
data=v2struct(load('behDataS.mat'));
%data.L1Error = lag(data.Error, 1);
data.L1Error(1) = nan; 
data.L1Error(2:length(data.L1Error)) = data.Error(1:(length(data.Error)-1)); 
%data=ConflictFlanker_ProcessDSVars('BEH_BC',data);
saveMatFile(data,'refinedDataS','refinedDataS');

%%

% % %Process2:Preparing all variables in dataset array
cd(BASIC.ALLGRANDIR)
data=v2struct(load('refinedDataS.mat'));
%data=ConflictFlanker_ProcessDSVars('CUTOFF',data);
data.MCOND=strcat(num2str(data.Target), "_", num2str(data.Flanker));
saveMatFile(data,'processedDataS','processedDataS');

%%
%AFTTER PREPROCESS!========================================================
%==========================================================================


%%
% % Process3:Aggregate data and get averages of interest!


%%
%FOR R !===================================================================
%==========================================================================
% acc=grpstats(data,{'SUBID'},{'mean'},'datavars','ACC')
% kvalue=grpstats(data,{'SUBID'},{'mean'},'datavars','mean_KVALUE')
% corr(kvalue.mean_mean_KVALUE,acc.mean_ACC);


%%
%PLOTTING!=================================================================
%==========================================================================


%%
%CORRELATION!==============================================================
%==========================================================================

%%


%%
% Done!
disp('Done!');

end