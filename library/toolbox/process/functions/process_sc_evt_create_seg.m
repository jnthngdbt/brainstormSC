function varargout = process_sc_evt_create_seg( varargin )
% process_sc_evt_create_seg: Create a new event group corresponding to a
% time segmentation with overlap
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'EVENTS > Create segmentation event group';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(435);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    % 
    sProcess.options.eventname.Comment = 'Event group name:';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'SEG';
    %
    sProcess.options.duration.Comment = 'Duration: ';
    sProcess.options.duration.Type    = 'value';
    sProcess.options.duration.Value   = {4,'sec',2};
    %
    sProcess.options.overlap.Comment = 'Overlap: ';
    sProcess.options.overlap.Type    = 'value';
    sProcess.options.overlap.Value   = {50,'%',0};
    %
    sProcess.options.replace.Comment = 'Replace existing events (same name)';
    sProcess.options.replace.Type    = 'checkbox';
    sProcess.options.replace.Value   = 1;
    % 
    sProcess.options.boxevt.Comment = 'Use those events as boxes in which to create segmentation(s): ';
    sProcess.options.boxevt.Type    = 'text';
    sProcess.options.boxevt.Value   = '';
    % 
    sProcess.options.label1.Comment = '<HTML>- If no box events are used, a unique segmentation covering the entire data is created.';
    sProcess.options.label1.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    
% Input options
eventname = sProcess.options.eventname.Value;    % New event group name
duration  = sProcess.options.duration.Value{1};  % Segments duration (s)
overlap   = sProcess.options.overlap.Value{1}/100; % Segments overlap (-)
isreplace = sProcess.options.replace.Value;      % Replaced same named evts
boxevt    = ExtractEvtNames(sProcess.options.boxevt.Value); 

% Initialize progress bar
if bst_progress('isVisible'), bst_progress('set', bst_progress('get')); end
bst_progress('start', 'Processing', 'Creating segmentation events...', 0, numel(sInputs));

for iInput=1:numel(sInputs)
    bst_progress('set', iInput);
    
    % Load data structure, time window and events informations
    DataMat = in_bst_data(sInputs(iInput).FileName);
    isRaw = strcmpi(sInputs(iInput).FileType, 'raw');
    if isRaw
        sfreq = DataMat.F.prop.sfreq;
        DataTimeWindow = [DataMat.F.prop.times(1); DataMat.F.prop.times(end)];
        sEvents = DataMat.F.events;
    else
        sfreq = 1/(DataMat.Time(2)-DataMat.Time(1));
        DataTimeWindow = [min(DataMat.Time); max(DataMat.Time)];
        sEvents = DataMat.Events;
    end
    
    % If no box events specified, the only timewindow is the entire data
    if isempty(boxevt)
        timewindow = DataTimeWindow;
    else
    % Otherwise, multiple timewindows for all specified box events
        DataEvtNames = {sEvents.label};
        timewindow = [];
        for iEvt=1:numel(boxevt)
            % Find wanted events
            idx = strcmpi(DataEvtNames,boxevt{iEvt});
            if ~any(idx)
                strMsg = sprintf('No events named [%s]',boxevt{iEvt});
                bst_report('Warning', sProcess, sInputs(iInput), strMsg);
                continue;
            end
            % Make sure it is extended
            EvtTimes = sEvents(idx).times;
            if size(EvtTimes,1)==1
                % If not, convert to extended by taking intervals
                strMsg = sprintf('Events [%s] are not extended. Creating extension as intervals.',boxevt{iEvt});
                bst_report('Warning', sProcess, sInputs(iInput), strMsg);
                EvtTimes = [EvtTimes; EvtTimes];
                EvtTimes(2,1:end-1) = EvtTimes(1,2:end);
                EvtTimes(2,end) = DataTimeWindow(end);
            end
            timewindow = [timewindow, EvtTimes];
        end
    end
    
    % Remove time windows smaller than wanted segments
    timewindow(:, diff(timewindow)<duration) = [];
    % Error if time windows are too small for wanted segments
    if isempty(timewindow)
        strMsg = 'Wanted segments duration are too big for available data';
        bst_report('Error', sProcess, sInputs(iInput), strMsg);
    end

    % Create segmentation times
    times = [];
    for iWin = 1:size(timewindow,2)
        WinTimes = timewindow(1,iWin):(1-overlap)*duration:(timewindow(2,iWin)-duration);
        WinTimes = [WinTimes;WinTimes+duration];
        % Add to main container
        times = [times, WinTimes];
    end
            
    % Create segmentation events
    newEvents = db_template('event');
    newEvents.label = eventname;
    newEvents.color = [0,0,0];
    newEvents.epochs = ones(1,size(times,2));
    newEvents.samples = round(times*sfreq);
    newEvents.times = times;
     
    % Report changes in .mat structure
    if isRaw
        if isreplace
            DataMat.F.events(strcmpi({DataMat.F.events.label},eventname)) = [];
        end
        DataMat.F.events = [DataMat.F.events,newEvents];
    else
        if isreplace
            DataMat.Events(strcmpi({DataMat.Events.label},eventname)) = [];
        end
        DataMat.Events = [DataMat.Events, newEvents];
    end
    % Save file definition
    bst_save(file_fullpath(sInputs(iInput).FileName), DataMat, 'v6',0);
end

OutputFiles = {sInputs.FileName};

end

% -------------------------------------------------------------------------
function EvtNames = ExtractEvtNames(EvtStr)

EvtNames = EvtStr;
if ~isempty(EvtNames)
    EvtNames  = textscan(EvtNames,'%s','delimiter',','); 
    EvtNames = EvtNames{1};
end
    
end
