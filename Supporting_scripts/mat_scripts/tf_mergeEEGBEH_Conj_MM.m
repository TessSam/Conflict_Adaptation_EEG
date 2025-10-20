function [DATA,EEG]=tf_mergeEEGBEH_Conj_MM(SUBID,DATA,FILEPATH,EXCON,DEPV)
% This merging function is assuming beh data and EEG data has same number
% of epoch. Check safer matching (but slow) for alternative ways.
% =========================================================================
keyboard;

tic
% Load EEG data for the subject
FILEPATH.SUBJECT=uniformLabel({'A','0',3},SUBID,'first');
FILEPATH.DATASETDIR=strcat(FILEPATH.EEGDIR,'/',FILEPATH.SUBJECT{:},'/DATASETS');
% FILEPATH.DATASETDIR=strcat(FILEPATH.EEGDIR,'/',FILEPATH.SUBJECT{:},'/DATASETS2');
cd(FILEPATH.DATASETDIR)%newer

% Load data and reshape it!
if EXCON.STIMLOCKED
    setFP=strcat(pwd,filesep,FILEPATH.SUBJECT,'_EEG_set.mat');
else
    setFP=strcat(pwd,filesep,FILEPATH.SUBJECT,'_R_EEG_set.mat');
end

%setFP=strcat(pwd,filesep,ls('*EEG_set.mat'));
for d=1:length(DEPV)
    if EXCON.STIMLOCKED
        dataFP(d)=strcat(pwd,filesep,FILEPATH.SUBJECT,'_EEG_set_',DEPV{d},'.bin');
    else
        dataFP(d)=strcat(pwd,filesep,FILEPATH.SUBJECT,'_R_EEG_set_',DEPV{d},'.bin');
    end
end

EEG=tf_loadbins(setFP{1},dataFP,DEPV);

% Change the field name depending on data type!
% switch DEPV{:},case{'eegraw'},fAdd='BINEPOCH';otherwise,fAdd='BINEPOCH.WAVELET';end

% Merge to Behavioral data
numepoch=EEG.BINEPOCH.numepochs;
evinfo=horzcat(EEG.EVENTLIST.eventinfo.code);
blcode=EXCON.BLOCKCODE;trcode=EXCON.TRIALCODE;
%sca;keyboard;
DATA=DATA(DATA.SUBID==SUBID,:);

% Prepare List of indexes
blInd=evinfo(ismember(evinfo,unique(DATA.Block)+blcode(1)))-blcode(1)';
trInd=evinfo(ismember(evinfo,unique(DATA.Trial)+trcode(1)))-trcode(1)';
epochInd=find(bsxfun(@(b,t)b&t,ismember(DATA.Block,blInd),ismember(DATA.Trial,trInd)));
control=false(1,3);%initialize to not to load anything,then change!
if any(strcmpi(DEPV,'eegraw')),eegraw=EEG.BINEPOCH.eegraw;control(1)=true;end
if any(strcmpi(DEPV,'eegpower')),eegpw=EEG.BINEPOCH.WAVELET.eegpower;control(2)=true;end
if any(strcmpi(DEPV,'eegphase')),eegph=EEG.BINEPOCH.WAVELET.eegphase;control(3)=true;end

rejectC=EEG.BINEPOCH.reject';
reject=nan(size(DATA,1),1);eegRaw=cell(size(DATA,1),1);
eegpwC=cell(size(DATA,1),1);eegphC=cell(size(DATA,1),1);

%keyboard;

switch size(DATA,1)==numepoch
    case {true},%liberal and faster way
        for be=1:numepoch
            beInd=epochInd(be);reject(beInd)= rejectC(be);
            if control(1),eegRaw(beInd)={eegraw(:,:,be)};end
            if control(2),eegpwC(beInd)={eegpw(:,:,:,be)};end
            if control(3),eegphC(beInd)={eegph(:,:,:,be)};end
        end
    case {false},%safe and slow way
        %keyboard;%**********
        disp('Epoch number mismatched! Check if its meant to be!')
        disp('Data will be merged based on event codes!!')

        %Get more information
        allevs=horzcat(EEG.EVENTLIST.eventinfo.bepoch);
        evinfo=horzcat(EEG.EVENTLIST.eventinfo.code);
        blcode=EXCON.BLOCKCODE;trcode=EXCON.TRIALCODE;
        numblock=DATA.Block;numtrial=DATA.Trial;
%keyboard;
        %Check event code and block/trial codes from beh data...
        for be=1:numepoch
            %keyboard;
            %Identify epoch to pull the data from!
            t_lock_ev=find(allevs==be);
            bl=evinfo(t_lock_ev+blcode(2));
            tr=evinfo(t_lock_ev+trcode(2));%beInd is computed from EVcodes
            beInd=find(numblock==bl-blcode(1)&numtrial==tr-trcode(1),1);
            
            %only trial is identified...
            if ~isempty(find(beInd,1)),
                
                reject(beInd)= rejectC(be);
                if control(1),eegRaw(beInd)={eegraw(:,:,be)};end
                if control(2),eegpwC(beInd)={eegpw(:,:,:,be)};end
                if control(3),eegphC(beInd)={eegph(:,:,:,be)};end
            end
            %keyboard;
        end
end
keyboard;
% Store everything Artifact rejection

DATA.eegraw=eegRaw;DATA.eegpower=eegpwC;DATA.eegphase=eegphC;DATA.reject=reject;
DATA=DATA(DATA.reject~=1,:);
%keyboard;
%DATA.reject(isnan(DATA.reject),:)=0
%DATA=DATA(~isnan(DATA.reject),:);%ask Melissa if this is necessary

%Erase unnecessary studd from EEG structure
EEG.EXCON=EXCON;
EEG.EVENTLIST=[];
EEG.BINEPOCH.eegraw=[];
EEG.BINEPOCH.WAVELET.eegpower=[];
EEG.BINEPOCH.WAVELET.eegphase=[];

%sca;keyboard;

% Catch some error
%sca;keyboard;
if size(DATA,1)==0,error('merging failed! Check merge event codes!');end
% if ~istable(DATA), DATA=dataset2table(DATA);end;%convert it to table!
toc

end
