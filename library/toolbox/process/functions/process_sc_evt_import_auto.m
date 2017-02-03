function varargout = process_sc_evt_import_auto( varargin )
% process_sc_evt_import_auto: Import a list of events files automaticaly
% using BST <process_evt_import>, finding the right raw input file in which
% to import based on each events file name (the target raw file name must 
% be part of the events file name)
% 
% Example:
% 
% Input raw files:
% - fileA.eeg
% - fileB.eeg
% - fileC.eeg
% 
% Selected event files:
% - fileB_spikes.vmrk
% - fileA_spikes.vmrk
% - fileC_spikes.vmrk
% 
% From the names, it is possible to match each event file to a raw file.
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'IMPORT > EVENTS > From file (name-matched raw-event files)';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(55);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    % File selection options
    SelectOptions = {...
        '', ...                               % Filename
        '', ...                               % FileFormat
        'open', ...                           % Dialog type: {open,save}
        'Import events...', ...               % Window title
        'ImportData', ...                     % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'multiple', ...                         % Selection mode: {single,multiple}
        'files', ...                          % Selection mode: {files,dirs,files_and_dirs}
        bst_get('FileFilters', 'events'), ... % Get all the available file formats
        'EventsIn'};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn
    % Option: MRI file
    sProcess.options.evtfile.Comment = 'Event file:';
    sProcess.options.evtfile.Type    = 'filename';
    sProcess.options.evtfile.Value   = SelectOptions;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    
    % Input options
    evtfiles = sProcess.options.evtfile.Value{1}; % Event files names
    format = sProcess.options.evtfile.Value{2};   % Event files format

    % Get raw files names
    for i=1:numel(sInputs)
        FileMat = in_bst_data(sInputs(i).FileName);

        % Raw filename parts
        [rawpath,rawname,rawext] = fileparts(FileMat.F.filename);
        if ~isempty(rawpath), rawpath = [rawpath,filesep]; end

        % Find match of list of input events files with current raw file 
        iMatch = find(cellfun(@(x)~isempty(strfind(upper(x),upper(rawname))),evtfiles));
        if isempty(iMatch)
            strMsg = 'No events file matched for this raw file';
            bst_report('Warning', sProcess, sInputs(i), strMsg)
            continue; 
        end
        
        % Import all matches
        for iEvt = 1:numel(iMatch)
            % Process: Events: Import from file
            sFiles = bst_process(...
                'CallProcess', 'process_evt_import', ...
                {sInputs(i).FileName}, [], ...
                'evtfile', {evtfiles{iMatch(iEvt)}, format});
        end
        
    end
    OutputFiles = {sInputs.FileName};
end



