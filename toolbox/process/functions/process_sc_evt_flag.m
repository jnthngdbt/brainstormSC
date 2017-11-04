function varargout = process_sc_evt_flag( varargin )
% process_sc_evt_flag: Flag an event group with another event group
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'EVENTS > Flag an event group with another event group';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(455);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw', 'data'};
    sProcess.OutputTypes = {'raw', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    % 
    sProcess.options.identified.Comment = 'Group A event group(s) name(s) (get flagged):';
    sProcess.options.identified.Type    = 'text';
    sProcess.options.identified.Value   = '';
    % 
    sProcess.options.identifier.Comment = 'Group B event group(s) name(s) (used to flag):';
    sProcess.options.identifier.Type    = 'text';
    sProcess.options.identifier.Value   = '';
    % -----------------------------------
    sProcess.options.sep1.Type = 'separator';
    sProcess.options.sep1.Comment = '';
    % 
    sProcess.options.label1.Comment = 'Flag method options: ';
    sProcess.options.label1.Type    = 'label';
    %
    sProcess.options.extension.Comment = {...
        '#1) Use group B extension and consider group A as punctual', ...
        '#2) Use group A extension and consider group B as punctual', ...
        '#3) Use both group A and B extensions (flag if overlap)'};
    sProcess.options.extension.Type    = 'radio';
    sProcess.options.extension.Value   = 1;
    % -----------------------------------
    sProcess.options.sep2.Type = 'separator';
    sProcess.options.sep2.Comment = '';
    % 
    sProcess.options.label2.Comment = 'If option #1 or #2 is selected: ';
    sProcess.options.label2.Type    = 'label';
    %
    sProcess.options.punctuality.Comment = {...
        'Convert extended events to punctual events using first sample', ...
        'Convert extended events to punctual events using middle sample', ...
        'Convert extended events to punctual events using last sample'};
    sProcess.options.punctuality.Type    = 'radio';
    sProcess.options.punctuality.Value   = 1;
% %     % -----------------------------------
% %     sProcess.options.sep3.Type = 'separator';
% %     sProcess.options.sep3.Comment = '';
% %     % 
% %     sProcess.options.label3.Comment = 'If option #3 is selected: ';
% %     sProcess.options.label3.Type    = 'label';
% %     %
% %     sProcess.options.overlap.Comment = {...
% %         'Flag if a group A event and B events overlap', ...
% %         'Convert extended events to punctual events using middle sample', ...
% %         'Convert extended events to punctual events using last sample'};
% %     sProcess.options.overlap.Type    = 'radio';
% %     sProcess.options.overlap.Value   = 1;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    
% INPUT OPTIONS

% Group A: Identified event group(s) name(s) (get flagged)
EvtAName = ExtractEvtNames(sProcess.options.identified.Value);
% Group B: Identifier event group(s) name(s) (used to flag)
EvtBName = ExtractEvtNames(sProcess.options.identifier.Value);
% Extension
switch sProcess.options.extension.Value
    case 1, extension = 'identifier';
    case 2, extension = 'identified';
    case 3, extension = 'both';
end
% Punctuality
switch sProcess.options.punctuality.Value
    case 1, punctuality = 'first';
    case 2, punctuality = 'middle';
    case 3, punctuality = 'last';
end

% Initialize progress bar
% % if bst_progress('isVisible'), bst_progress('set', bst_progress('get')); end
nProgress = numel(EvtAName)*numel(EvtBName)*numel(sInputs);
bst_progress('start', 'Processing', 'Flagging events...', 0, nProgress);

for iInput=1:numel(sInputs)
    
    
    % Load data structure and time window informations
    DataMat = in_bst_data(sInputs(iInput).FileName);
    isRaw = strcmpi(sInputs(iInput).FileType, 'raw');
                    
    % Compute for all identifier/identified pairs
    iProgress = 1;
    for iEvtA = 1:numel(EvtAName)
        for iEvtB = 1:numel(EvtBName)
            bst_progress('set', iProgress); iProgress = iProgress+1;
            if isRaw
                % Get both sets of events
                sEventA = DataMat.F.events(strcmpi({DataMat.F.events.label},EvtAName{iEvtA}));
                sEventB = DataMat.F.events(strcmpi({DataMat.F.events.label},EvtBName{iEvtB}));
                if isempty(sEventA) || isempty(sEventB), continue; end
                % Flag events
                newEvents = Compute(sEventA,sEventB,extension,punctuality);
                if isempty(newEvents), continue; end
                % Put new events in data file (overwriting pre-existing)
            	DataMat.F.events(strcmpi({DataMat.F.events.label},newEvents.label)) = [];
                DataMat.F.events = [DataMat.F.events,newEvents];
            else
                % Get both sets of events
                sEventA = DataMat.Events(strcmpi({DataMat.Events.label},EvtAName{iEvtA}));
                sEventB = DataMat.Events(strcmpi({DataMat.Events.label},EvtBName{iEvtB}));
                if isempty(sEventA) || isempty(sEventB), continue; end
                % Flag events
                newEvents = Compute(sEventA,sEventB,extension,punctuality);
                if isempty(newEvents), continue; end
                % Put new events in data file (overwriting pre-existing)
            	DataMat.Events(strcmpi({DataMat.Events.label},newEvents.label)) = [];
                DataMat.Events = [DataMat.Events, newEvents];
            end
            
        end
    end
    
    % Save file definition
    bst_save(file_fullpath(sInputs(iInput).FileName), DataMat, 'v6',0);
end

OutputFiles = {sInputs.FileName};

end

% -------------------------------------------------------------------------
function sEventC = Compute(sEventA,sEventB,extension,punctuality)
% Flag Brainstorm events sEventA (identified) with Brainstorm events 
% sEventB (identifier). 
% 

% Default inputs
if nargin<3, extension = 'both'; end
if nargin<4, punctuality = 'middle'; end
% Input checks
if numel(sEventA)>1 || numel(sEventB)>1
    error('Cannot process multiple events groups')
end

% Make sure the events are not empty, otherwise quit
sEventC = [];
if isempty(sEventA.samples) || isempty(sEventB.samples), return; end

% Make sure events have start and end points
StartA = sEventA.samples(1,:);
if size(sEventA.samples,1)==1, EndA = StartA;
else                           EndA = sEventA.samples(2,:);
end
StartB = sEventB.samples(1,:);
if size(sEventB.samples,1)==1, EndB = StartB;
else                           EndB = sEventB.samples(2,:);
end

% Extension and punctuality of events
switch lower(extension)
    case 'identified' 
    % --- Keep identified (A) extension and punctialize identifier (B)
        switch lower(punctuality)
            case 'first',  EndB = StartB;
            case 'middle', EndB = round(mean([StartB;EndB])); StartB = EndB;
            case 'last',   StartB = EndB;
        end
    % --- Keep identifier (B) extension and punctialize identified (A)
    case 'identifier'
        switch lower(punctuality)
            case 'first',  EndA = StartA;
            case 'middle', EndA = round(mean([StartA;EndA])); StartA = EndA;
            case 'last',   StartA = EndA;
        end
    % ---
    case 'both'
        % OK
end

% Create a flag vector
N = numel(StartA);
flag_i = zeros(1,N);
for i=1:N
    flag_i(i) = sum((StartA(i) <= EndB) & (EndA(i) >= StartB));
end

% Create flag-created events
sEventC = db_template('event');
sEventC.label = [sEventA.label,'_',sEventB.label];
sEventC.color = sEventA.color;
sEventC.samples = sEventA.samples(:,flag_i>0);
sEventC.times   = sEventA.times(:,flag_i>0);
sEventC.epochs  = sEventA.epochs(flag_i>0);

end

% -------------------------------------------------------------------------
function EvtNames = ExtractEvtNames(EvtStr)

EvtNames = EvtStr;
if ~isempty(EvtNames)
    EvtNames  = textscan(EvtNames,'%s','delimiter',','); 
    EvtNames = EvtNames{1};
end
    
end


