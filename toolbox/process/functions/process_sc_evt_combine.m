function varargout = process_sc_evt_combine( varargin )
% process_sc_evt_combine: Combine multiple event groups into another event
% group
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'EVENTS > Combine groups into another group';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(445);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    % 
    sProcess.options.eventname.Comment = 'Events to group:';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = '';
    % 
    sProcess.options.eventgroup.Comment = 'Combined group name: ';
    sProcess.options.eventgroup.Type    = 'text';
    sProcess.options.eventgroup.Value   = '';
    %
    sProcess.options.replace.Comment = 'Replace combined group (if already existing)';
    sProcess.options.replace.Type    = 'checkbox';
    sProcess.options.replace.Value   = 1;
%     % 
%     sProcess.options.label1.Comment = '<HTML>- If no box events are used, a unique segmentation covering the entire data is created.';
%     sProcess.options.label1.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    
% Input options
eventname  = ExtractEvtNames(sProcess.options.eventname.Value); 
eventgroup = sProcess.options.eventgroup.Value;

if ~sProcess.options.replace.Value
    % Add the combined group to the list of events to group, so that if it
    % already exists, new combined events will be added to it
    eventname = [eventname(:); {eventgroup}];
end

% Initialize progress bar
if bst_progress('isVisible'), bst_progress('set', bst_progress('get')); end
bst_progress('start', 'Processing', 'Creating combining events...', 0, numel(sInputs));

for iInput=1:numel(sInputs)
    bst_progress('set', iInput);
    
    % Load data structure, time window and events informations
    DataMat = in_bst_data(sInputs(iInput).FileName);
    isRaw = strcmpi(sInputs(iInput).FileType, 'raw');
    if isRaw
        sDataEvents = DataMat.F.events;
    else
        sDataEvents = DataMat.Events;
    end
    
    DataEvtNames = {sDataEvents.label}; % All event names in data
    isAllExtended = true;               % All events are extended
    sGrpEvents = [];                    % Grouped events container
    for iEvt=1:numel(eventname)
        % Find wanted events
        idx = strcmpi(DataEvtNames,eventname{iEvt});
        if ~any(idx) && ~strcmpi(eventname{iEvt},eventgroup)
            strMsg = sprintf('No events named [%s]',eventname{iEvt});
            bst_report('Warning', sProcess, sInputs(iInput), strMsg);
            continue;
        end
        % Check if extended
        isAllExtended = isAllExtended && size(sDataEvents(idx).times,1)==2;
        % Add new events info to container
        sGrpEventsNew.times   = sDataEvents(idx).times;
        sGrpEventsNew.samples = sDataEvents(idx).samples;
        sGrpEventsNew.epochs  = sDataEvents(idx).epochs;
        sGrpEvents = [sGrpEvents, sGrpEventsNew];        
    end
                
    % Create combined events group
    newEvents = db_template('event');
    newEvents.label   = eventgroup;
    newEvents.color   = [0,0,0];
    if ~isempty(sGrpEvents)
        newEvents.epochs  = [sGrpEvents.epochs];
        if isAllExtended
            newEvents.samples = [sGrpEvents.samples];
            newEvents.times   = [sGrpEvents.times];
        else
            newEvents.samples = [sGrpEvents.samples(1,:)];
            newEvents.times   = [sGrpEvents.times(1,:)];
        end
    else
        strMsg = sprintf('No events were combined');
        bst_report('Warning', sProcess, sInputs(iInput), strMsg);
    end
    
    % Sort by position
    [tmp,iSort] = sort(newEvents.times(1,:),'ascend');
    newEvents.times   = newEvents.times(:,iSort);
    newEvents.samples = newEvents.samples(:,iSort);
    newEvents.epochs  = newEvents.epochs(:,iSort);
    
    % Report changes in .mat structure
    if isRaw
        DataMat.F.events(strcmpi({DataMat.F.events.label},newEvents.label)) = [];
        DataMat.F.events = [DataMat.F.events,newEvents];
    else
        DataMat.Events(strcmpi({DataMat.Events.label},newEvents.label)) = [];
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
