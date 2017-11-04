function varargout = process_sc_fasst_score_import( varargin )
% process_sc_fasst_score_import: Return sleep scoring done in FASST.
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'SLEEP > FASST > Import Sleep Scoring (unstable)';
    sProcess.FileTag     = '';
    sProcess.Description = '<HTML>Import sleep scoring in an existing FASST file';
    sProcess.Category    = 'custom';
    sProcess.Index       = sc_bst_process_index(2144);
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
    
    bst_progress('start', 'FASST Sleep Scoring', 'Importing');
    
    for iFile = 1:length(sInputs)

        % ===== GET DATA =====
        % Load the raw file descriptor
        DataMat = in_bst_data(sInputs(iFile).FileName);
        sFile = DataMat.F;
        
        % ===== IMPORT SLEEP SCORING FROM FASST FILE =====
        [p,n,e] = fileparts(sFile.filename); p = [p,filesep];
        fasstFile = [p,n,'.mat'];
        if exist(fasstFile,'file')==0
            error('File must be scored using process [FASST: Score/Review Sleep]');
        end
        % Create BST importable file
        EventFile = export_scoring_bst({fasstFile});
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

% =========================================================================
function varargout = export_scoring_bst(filenames)

% BRAINSTORM STRUCT (array of struct for different labels)
%          label: 'SPIN'
%          color: [0 1 1]
%         epochs: [1 1 1 1 1]
%        samples: [2x5 double]
%          times: [2x5 double]
%     reactTimes: []
%         select: 1

% Get filenames
if nargin==0, filenames = get_files; end
if ~iscell(filenames), filenames = {filenames}; end

% Colors are extracted from FASST
scoreLabel = {'WAKE','NREM1','NREM2','NREM3','NREM4','REM','MVT','NA'};
scoreColor = [[0.2 0.75 0.6]; [0 0.8 1]; [0.1 0.5 0.9]; [0.1 0.2 0.8]; ...
              [0.1 0.15 0.5]; [0.5 0.5 0.9]; [0.9 0.4 0.4]; [0.9 0.6 0.3]];

Nf = numel(filenames);
status.isScored = ones(1,Nf);
for i=1:Nf
    % Load the file
    [p n e] = fileparts(filenames{i}); p = [p filesep];
    load([p,n,e]); % D
    
    if ~isfield(D.other,'CRC'), status.isScored(i) = 0; continue; end
    if ~isfield(D.other.CRC,'score'), status.isScored(i) = 0; continue; end

    
    fs  = D.Fsample;                % Sampling frequency
    dt  = 1/fs;                     % Sampling time increment
    scr = D.other.CRC.score{1,1};   % Sleep scores vector
    win = D.other.CRC.score{3,1};   % Sleep scoring window (s)
    win_i = round(win*fs);          % Sleep scoring window (samp)
    Ns  = numel(scr);               % Number of sleep scores
    pos_i = (0:Ns-1)*win_i + 1;     % Positions in file (samp)
    pos = pos_i/fs;                 % Positions in file (s)
    
    scr(isnan(scr)) = 7; % NA
    
    samples = [pos_i; pos_i+win_i-1];
    times = [pos; pos+win-dt];
    
    % Create the BST array of events structs
    events = []; % Will contain the BST events structs
    uScr = unique(scr); 
    for k=1:numel(uScr)
        
        iScr = scr==uScr(k);
        Ni = sum(iScr);
        
        events_i.label      = scoreLabel{uScr(k)+1};
        events_i.color      = scoreColor(uScr(k)+1,:);
        events_i.epochs     = ones(1,Ni);
        events_i.samples    = samples(:,iScr);
        events_i.times      = times(:,iScr);
        events_i.reactTimes = [];
        events_i.select     = 1;
        
        events = [events, events_i];
    end
    
    outFiles{i} = [p,'events_fasst_scoring_',n,'.mat'];
    save(outFiles{i},'events');
end

if nargout>=1, varargout{1} = outFiles; end

end

function filenames = get_files

[filenames, pname] = uigetfile( ...
{'*.mat' , 'MAT Files (*.mat)'}, ...
'Pick FASST MAT files', ...
'MultiSelect', 'on');

if ~iscell(filenames), filenames = {filenames}; end

N = numel(filenames);
for i=1:N
    filenames{i} = [pname,filenames{i}];
end

end






