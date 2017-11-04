function varargout = process_sc_fasst_score( varargin )
% process_sc_fasst_score: Launch FASST toolbox for sleep scoring. 

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'SLEEP > FASST > Score/Review Sleep (unstable)';
    sProcess.FileTag     = ' | scoring';
    sProcess.Description = '<HTML>Launch FASST for sleep scoring';
    sProcess.Category    = 'custom';
    sProcess.Index       = sc_bst_process_index(2142);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    
    bst_progress('start', 'FASST Sleep Scoring', 'Preparing');
    
    for iFile = 1:length(sInputs)

        % ===== GET DATA =====
        % Load the raw file descriptor
        DataMat = in_bst_data(sInputs(iFile).FileName);
        sFile = DataMat.F;
        sChan = in_bst_channel(sInputs(iFile).ChannelFile);
        
        % ===== LAUCH FASST FOR SLEEP SCORING =====
        if isempty(which('crc_main'))
            error('FASST must be installed');
        end
        [p,n,e] = fileparts(sFile.filename); p = [p,filesep];
        switch sFile.device 
            case 'BRAINAMP'
                fasstFile = [p,n,'.vhdr'];
            otherwise
                fasstFile = [p,n,e];
        end
        % Check extension for EDF
        switch lower(e)
            case '.rec'
                error('You must rename file extensions (.rec -> .edf)')
        end
        
        % Following is extracted from DIS_SELCHAN (95)
        D = crc_eeg_load(fasstFile);
        i = 1;
 % %         handles.chan{i} = upper(chanlabels(D));
        chan = {'F3','Fz','F4','C3','Cz','C4','P3','Pz','P4','O1','Oz','O2'};
        handles.chan{i} = chan;
        handles.index = [];
        for k=1:numel(chan)
            handles.index = [handles.index,find(strcmpi(chan{k},chanlabels(D)))];
        end
        handles.file{i} = fullfile(D.path,D.fname);
        handles.Dmeg{i} = D;
        handles.date{i} = zeros(1,2);
        handles.date{i}(1) = datenum([D.info.date D.info.hour]);
        handles.date{i}(2) = handles.date{i}(1) + ...
        datenum([ 0 0 0 crc_time_converts(nsamples(D)/fsample(D))] );
        handles.dates(i,:) = handles.date{i}(:);        
        % Call for channel selection; from there, select SCORE
        dis_selchan(handles);
        
    end
    % Return all the input files
    OutputFiles = {sInputs(iFile).FileName};
            
    bst_progress('stop');
end




