function [bData,eData,headerDim]=tf_trans2env_MM(DATA,EEG,FILEPATH,S)
% save trial-by-trial EEG data intp h5file to be transfered to other
% enviornment such as R and Python.

% -OPTIONS:
% depV= {'eegraw','eegpower','eegphase'};
% fband = {whatever frequency band,'ALL','NONE'};
% elec = {whatever range of electrodes to include, 'ALL'}
% zscore = apply z-scoreing among electrodes or not


% =========================================================================

% Get basic info and reshape EEG
timeL=EEG.BINEPOCH.times;
freqL=EEG.BINEPOCH.WAVELET.frex;
chanL={EEG.BINEPOCH.chanlocs.labels};%all electrodes
chanLH=chanL(~ismember(chanL,{'HEOG','VEOG'}));%whead electrodes
%keyboard;
eeg=unpackcell(DATA.(S.depV{1}));%elec,freq,time,trial
numD=ndims(eeg);

% Adjust dimention of the data (trial,time,fb,elec)
i=[];if strcmpi(S.depV{:},'eegraw'),i=4;end
insert=@(a,x,n)cat(2,x(1:n),a,x(n+1:end));
sv=insert(i,fliplr(1:numD),3);eeg=permute(eeg,[sv]);
sprintf('The original dimension is %d trial, %d  time points, %d fb, %d electrodes!!',size(eeg))

% Select time and electrodes
if ~isfield(S,'timeC'),timeC=minmax(timeL);end;timeV=dsearchn(timeL',timeC')';
if strcmpi(S.elec,'ALL'),chanV=1:length(chanLH);else chanV=find(strcmpi(chanLH,S.elec));end
%keyboard;
%eeg=idxDR(eeg,[2,4],{timeV(1):timeV(2),chanV});
eeg=idxDR(eeg,[2,4],{timeV(1):timeV(2),chanV});

% Average across frequency band (for eegpower and eegphase)
if ~all(strcmpi(S.fband,'NONE'))
    for fb=1:length(S.fband)
        %Get the index of frequency values
        fbN=S.fband{fb};freqT=EEG.FREQINFO.(fbN);
        fr=dsearchn(freqL',minmax(freqT)')';
        d=mean(eeg(:,:,fr(1):fr(2),:),3);
        if S.zscore,d=zscore(d,0,4);end
        eegAV(:,:,fb,:)=d;
    end
end

% Save as h5 file!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dpath=[FILEPATH.SUBDIR{S.subInd},filesep,'DATASETS'];
dfile=strcat(dpath,filesep,FILEPATH.SUBJECT{:},'_',S.fext{:},'.h5');
if any(exist(dfile,'file')),delete(dfile);end;cd(dpath);

%Write EEG data (trial,time,fb[this could be size of 1],elec)
h5create(dfile,[filesep,S.depV{:}],size(eegAV),'Datatype','single');
h5write(dfile,[filesep,S.depV{:}], eegAV);%h5disp(dfile)
sprintf('Saving data as %d trial, %d  time points, %d fb, %d electrodes!!',size(eegAV));

%Write BLOCK & TRIAL index
btrial=single(horzcat(DATA.Block,DATA.Trial));
% btrial=single(horzcat(DATA.BLOCK,DATA.TRIAL,DATA.EVNUM));
h5create(dfile,[filesep,'IDX'],size(btrial),'Datatype','single');
h5write(dfile,[filesep,'IDX'], btrial);%h5disp(dfile)

% Write dimension of each ??
% h5create(dfile,[filesep,'TrElFbTi'],size(btrial),'Datatype','single');

%Write header?? or other information??

end