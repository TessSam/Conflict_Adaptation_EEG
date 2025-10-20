function [EEG] = TimeFreq_EEGLAB_Preprocess_ConflictAdaptation(PREP)

%%

%convenient wrapper for doing preprocessing and artifact rejection of data
%from our lab.
% EEG = [];
% ALLEEG = [];
% CURRENTSET = [];
% params = [];
% p = inputParser;

%%
% -------------------------------------------------------------------------
% STEP1: Loading, Channel Setting, Calibration, Filtering
% -------------------------------------------------------------------------
% Create DATASETS folder Fire up EEGlab
cd(PREP.SUBDIR{PREP.i});
if ~exist(fullfile(cd, 'DATASETS'),'dir'),mkdir('DATASETS');end
PREP.DATASETDIR=strcat(PREP.SUBDIR{PREP.i},'/DATASETS');
setfiles=dir('*.vhdr');
PREP.SETFILES={setfiles.name};
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
%cd(PREP.SUBDIR);
%PREP.SETFILES={'Flank005.vhdr'};

%[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%% Loop over n SETFILES 
% When data are stored in separate files (e.g., subjects went to restroom etc), run preprocess independently.
% Each file has its own time samples and indexes so epoching must reflect differences. 

for si = 1:length(PREP.SETFILES)

    % %Import the raw data and make an eeglab dataset
    %setName = erase(PREP.SETFILES{si},'.vhdr');
    %EEG = pop_loadbv(PREP.SUBDIR, [PREP.SETFILES{si}]);
    EEG = pop_loadbv(PREP.SUBDIR{PREP.i}, [PREP.SETFILES{si}]);
    [ALLEEG EEG CURRENTSET]=pop_newset(ALLEEG,EEG,1,'gui','off','setname','raw_data');

    
    %%
    % %Define channel locations and information
    chTemp = loadbvef(PREP.ELECDIR);
    %chTemp = [chTemp,structfun(@(x)[],chTemp(end),'Uni',false)];% copy template for EOG
    %chTemp(end).labels = 'EOG';
    %chRemap = vertcat(PREP.ELECREMAP{:});
    %chRIDX = ~ismember({EEG.chanlocs.labels},{chTemp.labels});
    % chRIDX(strcmpi('EOG',{EEG.chanlocs.labels})) = false;% ignore EOG
    % NOTE: does it make sence to keep A1/A2 locations....???
    
%     for ch = 1:length({EEG.chanlocs.labels})
%         chL = EEG.chanlocs(ch).labels;
%         if chRIDX(ch),chL = chRemap{ismember({chRemap{:,2}},chL),1};end
%         ch_o = logical(ismember({chTemp.labels},chL));
%         ch_info = chTemp(ch_o);ch_info.labels = EEG.chanlocs(ch).labels;
%         EEG.chanlocs(ch) = catstruct(EEG.chanlocs(ch),ch_info);
%     end
    
    % %Define channel locations and information (old)
    % EEG=pop_chanedit(EEG,'load',{PREP.ELECDIR, 'filetype', 'autodetect'});
    % [ALLEEG EEG] = eeg_store(ALLEEG,EEG,CURRENTSET);
            
    %%
    % % Do average re-referecing??:
    % % convert common reference EEG data to average references (check bva_tools)
    %chIDX = find(ismember({EEG.chanlocs.labels},{'EOG','VEOG','HEOG'}));
    %if PREP.DOREREF,EEG.data = reref(EEG.data,[],'exclude',chIDX);end
    
    %%
    
    % % Do reference adjustment (for left-right mastoid REFs)!
    %  Check Chapter 3 (p44) equation for average mastoid REFs
    if PREP.DOREREF
        %refE1 = ismember({EEG.chanlocs.labels},'A1'); % reference will be omited
        %refE2 = ismember({EEG.chanlocs.labels},'A2');
        refE1 = ismember({EEG.chanlocs.labels},'TP9'); % reference will be omited
        refE2 = ismember({EEG.chanlocs.labels},'TP10');
        
        if (~isempty(find(refE2, 1)))
            data_refE2 = repmat(EEG.data(1,:),length(EEG.chanlocs),1);
            EEG.data = EEG.data - (1/2 * data_refE2); % a' = a - (r/2)
        end
    end
    
    % %%
    % % Do reference adjustment (old, specific to Oregon)!
    % if PREP.DOCALIB,[EEG]=eegm_calib(EEG,PREP.BOOTH);end
    
    %%
    % Filtering (make sure to do this to remove drifts)
    if PREP.DOFILTER
        EEG  = pop_basicfilter(EEG, PREP.FILTERCHAN , 'Boundary', 'boundary', 'Cutoff', PREP.FILTERCUT, ...
            'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
    end
    
    % Save the result into .set??
    if PREP.SAVESET,EEG=pop_saveset(EEG,'filename',[PREP.SUBLABEL{PREP.i} '_filtered.set'],'filepath',PREP.DATASETDIR);end
    
    %%
    % Filtering for EOG...(why does it have super high-freq noise?)
    chIDX = find(ismember({EEG.chanlocs.labels},{'Fp1','Fp2'}));
    %chIDX = find(ismember({EEG.chanlocs.labels},{'EOG','VEOG','HEOG'}));
    EEG  = pop_basicfilter(EEG,chIDX, 'Boundary', 'boundary', 'Cutoff',[0,8], ...
        'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
    
    
    %%
    % Re(down)sampling
    if PREP.DORESAMPLE,[EEG] = pop_resample(EEG, PREP.RESAMPLE);end
    
    %%
    % -------------------------------------------------------------------------
    % STEP2: Binning and Epoching
    % -------------------------------------------------------------------------
    % Create eventlist to prep for binning and epoching
    EEG=pop_creabasiceventlist(EEG,'AlphanumericCleaning','on','BoundaryNumeric',{-99},...
        'BoundaryString',{'b'},'EVENTLIST',[PREP.DATASETDIR filesep sprintf('%s_bin_eventlist.txt',PREP.SUBLABEL{PREP.i})] ); % GUI: 29-Apr-2016 13:50:13
    
    %%
    % Binning!
    if PREP.DOBIN,EEG=pop_binlister(EEG,'BDF',PREP.BDFDIR,'IndexEL',1,'SendEL2','EEG','Voutput','EEG');end;
    [ALLEEG EEG CURRENTSET]=pop_newset(ALLEEG,EEG,CURRENTSET,'setname','binned_data');
    
    % Save the result into .set??
    if PREP.SAVESET,EEG=pop_saveset(EEG,'filename',[PREP.SUBLABEL{PREP.i} '_binned.set'],'filepath',PREP.DATASETDIR);end
    
    %%
    % Epoching!
    if PREP.DOEPOCH,EEG=pop_epochbin(EEG,[PREP.EPOCHDUR{:}],PREP.BASELINE);end
    [ALLEEG EEG CURRENTSET]=pop_newset(ALLEEG,EEG,CURRENTSET,'setname','epoched_data');
    
    % Save the result into .set??
    if PREP.SAVESET,EEG=pop_saveset(EEG,'filename',[PREP.SUBLABEL{PREP.i} 'epoched.set'],'filepath',PREP.DATASETDIR);end
    
    %%
    % -------------------------------------------------------------------------
    % STEP3: Artifact rejection
    % -------------------------------------------------------------------------
    
    % Artifact rejection:Blocking
    if PREP.DOBLOCK,EEG=pop_artflatline(EEG,'Twindow',[PREP.AFEPOCHDUR{:}],'Threshold',PREP.FLATCUT,'Duration',PREP.FLATWIND,'Channel',PREP.FLATCHAN,'Flag',PREP.FLATFLAG);end
    
    % Artifact rejection:Eye blink
    if PREP.DOEYEBLINK,EEG=pop_artmwppth(EEG,'Twindow',[PREP.AFEPOCHDUR{:}],'Threshold',PREP.BLINKCUT,'Windowsize',PREP.BLINKWIND,'Windowstep',PREP.BLINKSTEP,'Channel',PREP.BLINKCHAN,'Flag',PREP.BLINKFLAG);end
    
    % Artifact rejection:Eye movement
    %if PREP.DOEYEMOVE,EEG=pop_artstep(EEG,'Twindow',[PREP.AFEPOCHDUR{:}],'Threshold',PREP.EYEMCUT,'Windowsize',PREP.EYEMWIND,'Windowstep',PREP.EYEMSTEP,'Channel',PREP.EYEMCHAN,'Flag',PREP.EYEMFLAG);end
    if PREP.DOEYEMOVE,EEG=pop_artmwppth(EEG,'Twindow',[PREP.AFEPOCHDUR{:}],'Threshold',PREP.EYEMCUT,'Windowsize',PREP.BLINKWIND,'Windowstep',PREP.BLINKSTEP,'Channel',PREP.EYEMCHAN,'Flag',PREP.EYEMFLAG);end
    
    %write out artifact rejection output
    pop_summary_AR_eeg_detection(EEG,[PREP.DATASETDIR filesep sprintf('%s_AR_Summary.txt',PREP.SUBLABEL{PREP.i})]);
    
    %%
    % Keep original results
    EEGKEEP(si) = {EEG};
    
end
% -------------------------------------------------------------------------
% STEP4: Save results!
% -------------------------------------------------------------------------

% Were there separated files of EEG?
if length(EEGKEEP)==1,EEG = EEGKEEP{:};else, EEG = tf_appendEEG_BP(EEGKEEP);end

% Save the result into .set??
if PREP.SAVESET,EEG=pop_saveset(EEG,'filename',[PREP.SUBLABEL{PREP.i} 'preprocessed.set'],'filepath',PREP.DATASETDIR);end

% Save most udated EEG structure as a mat file
cd(PREP.DATASETDIR);
saveMatFile(EEG,'EEG',strcat(PREP.SUBLABEL{PREP.i},'_','EEG_preprocessed'));
%saveMatFile(EEG,'EEG',strcat('A105_','EEG_preprocessed'));

% Give message
fprintf('~~~~~Finished pre-processing for %s !!~~~~~\n.',PREP.SUBLABEL{PREP.i});
%fprintf('~~~~~Finished pre-processing for %s !!~~~~~\n.',setName);

cd('../..')
end

% % Bind separated raw EEG files
% function EEGm = tf_appendEEG_BP(EEG)
% 
% % Master structure
% EEGm = EEG{1};
% 
% % Merge data...currently ignoring most of stuff in "reject" field
% for d = 2:length(EEG)
%     EEGm.trials = EEGm.trials + EEG{d}.trials;
%     EEGm.data = cat(3,EEGm.data,EEG{d}.data);
%     EEGm.epoch = [EEGm.epoch, EEG{d}.epoch];
%     EEGm.event = [EEGm.event, EEG{d}.event];
%     EEGm.urevent = [EEGm.urevent, EEG{d}.urevent];
%     EEGm.EVENTLIST.eventinfo = [EEGm.EVENTLIST.eventinfo, EEG{d}.EVENTLIST.eventinfo];
%     EEGm.EVENTLIST.trialsperbin = EEGm.EVENTLIST.trialsperbin + EEG{d}.EVENTLIST.trialsperbin;
%     EEGm.reject.rejmanual = [EEGm.reject.rejmanual, EEG{d}.reject.rejmanual];
%     EEGm.reject.rejmanualE = [EEGm.reject.rejmanualE, EEG{d}.reject.rejmanualE];
% end
% 
% end