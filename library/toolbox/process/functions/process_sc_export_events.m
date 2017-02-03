function varargout = process_sc_export_events( varargin )
% process_sc_export_events: Just a scripting way to export events using
% builtin function EXPORT_EVENTS. Exports in Brainstorm format.
% 
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'EXPORT > Events';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1755);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    % Specific event names to export
    sProcess.options.events.Comment = 'Names of specific events to export (empty=All): ';
    sProcess.options.events.Type    = 'text';
    sProcess.options.events.Value   = '';
    % Export format
    sProcess.options.format.Comment = 'Events export format: ';
    sProcess.options.format.Type    = 'combobox';
    sProcess.options.format.Value   = {1,{'BST (*.mat)','CTF (*.mrk)','FIF (*.eve)','ARRAY-TIMES (*.txt)'}};
    % Export filename suffix
    sProcess.options.suffix.Comment = 'Suffix to add to raw filename (ex: filename.[suffix].*): ';
    sProcess.options.suffix.Type    = 'text';
    sProcess.options.suffix.Value   = 'EVENTS';
    %
    sProcess.options.exportdir.Comment = {'Export in protocol directory','Export where each raw file is'};
    sProcess.options.exportdir.Type    = 'radio';
    sProcess.options.exportdir.Value   = 1;
%     % Path where to export
%     sProcess.options.exportdir.Comment = 'Events are exported in raw files directory.';
%     sProcess.options.exportdir.Type    = 'label';
%     sProcess.options.exportdir.Value   = '';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {sInputs.FileName}; % Return input file names

% Get current progressbar position
if bst_progress('isVisible')
    curProgress = bst_progress('get');
else
    curProgress = [];
end

% Suffix to identify exported events files
suffix = sProcess.options.suffix.Value;
if isempty(suffix), suffix = 'EVENTS'; end

% Specific events to export
evtNames = sProcess.options.events.Value;
if ~isempty(evtNames)
    evtNames = textscan(sProcess.options.events.Value,'%s','delimiter',',');
    evtNames = evtNames{1};
end

% Protocol directory
sProtocol = bst_get('ProtocolInfo');
protocolDir = [sProtocol.STUDIES,filesep];

% Loop through files
bst_progress('text', 'Exporting events...');
Ni = numel(sInputs);
for iInput=1:Ni
    bst_progress('set',  round(curProgress + 100*iInput/Ni));
        
    % Get raw file struct descriptor
    sMat = in_bst_data(sInputs(iInput).FileName);
    sFile = sMat.F;
    
    % Raw filename parts
    [rawPath,rawName,rawExt] = fileparts(sFile.filename);
    if ~isempty(rawPath), rawPath = [rawPath,filesep]; end

    % Determine export directory
    switch sProcess.options.exportdir.Value
        case 1 % Protocol directory
            exportDir = protocolDir;
        case 2 % Raw file directory
            exportDir = rawPath;
    end
    
    % Determine export format
    switch sProcess.options.format.Value{1}
        case 1 % '.mat',   FileFormat = 'BST';
            exportExt = '.mat';
        case 2 % '.mrk',   FileFormat = 'CTF';
            exportExt = '.mrk';
        case 3 % '.eve',   FileFormat = 'FIF';
            exportExt = '.eve';
        case 4 % '.txt',   FileFormat = 'ARRAY-TIMES';
            exportExt = '.txt';
    end
    
    % Create export file name
    exportFilename = [exportDir,'events.',rawName,'.',suffix, exportExt];
    
    % Only keep events to export and error checking
    if isempty(sFile.events)
        strMsg = 'No events in file!';
        bst_report('Warning', sProcess, sInputs(iInput), strMsg)
    end
    if ~isempty(evtNames)
        fileEvtNames = {sFile.events.label}; % All file events
        iEvt = ismember(fileEvtNames,evtNames); % Find specified event
        if ~any(iEvt)
            strMsg = 'Specified events are not present in file!';
            bst_report('Warning', sProcess, sInputs(iInput), strMsg)
        end
        % Remove non-selected events (it will not be saved in raw file!!!!)
        sFile.events(~iEvt) = [];
    end
    % Export
    export_events( sFile, exportFilename );
    
end

end
