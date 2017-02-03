function varargout = process_sc_fasst_spindle_import( varargin )
% 

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'SLEEP > FASST > Import Spindle Detections (unstable)';
    sProcess.FileTag     = '';
    sProcess.Description = '<HTML>Import spindles detections from an existing FASST file';
    sProcess.Category    = 'custom';
    sProcess.Index       = sc_bst_process_index(2148);
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
    
    bst_progress('start', 'FASST Spindle Detections', 'Importing');
    
    for iFile = 1:length(sInputs)

        % ===== GET DATA =====
        % Load the raw file descriptor
        DataMat = in_bst_data(sInputs(iFile).FileName);
        sFile = DataMat.F;
        
        % ===== IMPORT SLEEP SCORING FROM FASST FILE =====
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
            error('File must be scored and spindle detected in FASST');
        end
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




