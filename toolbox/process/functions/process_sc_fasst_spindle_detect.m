function varargout = process_sc_fasst_spindle_detect( varargin )
% 

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'SLEEP > FASST > Detect Spindles (unstable)';
    sProcess.FileTag     = '';
    sProcess.Description = '<HTML>Launch FASST automatic spindles detection from an existing FASST file';
    sProcess.Category    = 'custom';
    sProcess.Index       = sc_bst_process_index(2146);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    
    sProcess.options.bandfc.Comment = 'Frequency band: ';
    sProcess.options.bandfc.Type    = 'range';
    sProcess.options.bandfc.Value   = {[8,20], 'Hz',1};

    sProcess.options.scorer.Comment = 'Scorer ID number: ';
    sProcess.options.scorer.Type    = 'value';
    sProcess.options.scorer.Value   = {1,'',0};

    sProcess.options.stagesp.Comment = 'Consider only stage: ';
    sProcess.options.stagesp.Type    = 'value';
    sProcess.options.stagesp.Value   = {[2 3 4],'',0}; % ,'e.g.: [2,3] for NREM2 and NREM3'

    sProcess.options.wav.Comment = 'Antero-posterior classification';
    sProcess.options.wav.Type    = 'checkbox';
    sProcess.options.wav.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    
    bst_progress('start', 'FASST Spindle Detection', 'Detecting');
    
    for iFile = 1:length(sInputs)

        % ===== GET DATA =====
        % Load the raw file descriptor
        DataMat = in_bst_data(sInputs(iFile).FileName);
        sFile = DataMat.F;
        
        % ===== LAUNCH FASST SPINDLE DETECTION =====
        [p,n,e] = fileparts(sFile.filename); p = [p,filesep];
%         switch sFile.device 
%             case 'EGI'
%                 lm1020File = [p,n,'_LM1020','.mat'];
%                 fasstFile = lm1020File;
%             otherwise
%                 fasstFile = [p,n,'.mat'];
%         end
        fasstFile = [p,n,'.mat'];
        if exist(fasstFile,'file')==0
            error('File must be scored in FASST');
        end
        
        handles.highfc  = sProcess.options.bandfc.Value{1}(1);
        handles.lowfc   = sProcess.options.bandfc.Value{1}(2);
        handles.review  = 0;
        handles.fname   = fasstFile;
        handles.reref   = 1;
        handles.scorer  = sProcess.options.scorer.Value{1};
        handles.stagesp = sProcess.options.stagesp.Value{1};
        handles.analyse = 3; % Weird; not stage for threshold, its a type of analysis..
        handles.Begpts  = 1;
        handles.Endpts  = numel(DataMat.Time);
        handles.wav     = sProcess.options.wav.Value; % Wavelet Analysis (also ANT + POST)
        crc_SP_detect(handles)
        
        % Create BST importable file
        EventFile = fasst.sleep.export_spindle_bst('files',{fasstFile});
        % Import to BST
        [sFile, newEvents] = import_events(sFile, EventFile{1}, 'BST');
        % Only save changes if something was change
        if ~isempty(newEvents)
            DataMat.F = sFile;
            % Keep only time window
            DataMat.Time = DataMat.Time([1, end]);
            % Save file definition
            save(file_fullpath(sInputs(iFile).FileName), '-struct', 'DataMat');
            % Report number of detected events
            bst_report('Info', sProcess, sInputs(iFile), sprintf('Added to file: %d events in %d different categories', size([newEvents.times],2), length(newEvents)));
        else
            bst_report('Error', sProcess, sInputs(iFile), 'No events read from file.');
        end
        % Return all the input files
        OutputFiles = {sInputs(iFile).FileName};
    end
    
            
    bst_progress('stop');
end




