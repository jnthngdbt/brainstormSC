function varargout = process_sc_report_has( varargin )
% process_sc_report_has: Reports which datasets (raw files) have (or not) a
% specific feature (channel, event, ...).
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'REPORT > Database Has analysis';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(255);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    % Option: Has (or not) ALL or has (or not) ANY
    sProcess.options.rule.Comment = {'Report files that has (or not) ALL of what is specified','Report files that has (or not) ANY of what is specified'};
    sProcess.options.rule.Type    = 'radio';
    sProcess.options.rule.Value   = 1;
    %
    sProcess.options.chn.Comment = 'Channel name: ';
    sProcess.options.chn.Type    = 'text';
    sProcess.options.chn.Value   = '';
    %
    sProcess.options.evt.Comment = 'Event name: ';
    sProcess.options.evt.Type    = 'text';
    sProcess.options.evt.Value   = '';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

% === INPUT PARAMETERS ===

% Get channel names specified
in.channel = [];
if ~isempty(sProcess.options.chn.Value)
    in.channel = textscan(sProcess.options.chn.Value,'%s','delimiter',',');
    in.channel = in.channel{1};
end

% Get event names specified
in.event = [];
if ~isempty(sProcess.options.evt.Value)
    in.event = textscan(sProcess.options.evt.Value,'%s','delimiter',',');
    in.event = in.event{1};
end

% Rule function 
switch sProcess.options.rule.Value
    case 1, ruleFcn = @all;
    case 2, ruleFcn = @any;
end

% === COMPUTE ===

[has] = Has(sInputs, in, ruleFcn);

% === REPORT ===

Report(sInputs, sProcess, 'channel', in.channel, has.channel);
Report(sInputs, sProcess, 'event', in.event, has.event);


% === EXIT ===
OutputFiles = {sInputs.FileName}; % Return input file names

end

% -------------------------------------------------------------------------
function [has] = Has(sInputs, in, ruleFcn)
% Brainstorm must be running

if nargin<3, ruleFcn = @all; end

Ni = numel(sInputs);
has.channel = false(1,Ni);
has.event = false(1,Ni);
bst_progress('start', 'Processing', 'Extracting "Has" information from files...', 0, Ni);
for iInput=1:Ni
    bst_progress('set', iInput);
    
    % Get data and channel structs
    DataMat = in_bst_data(sInputs(iInput).FileName,sInputs(iInput).ChannelFile);
    ChanMat = in_bst_channel(sInputs(iInput).ChannelFile);
         
    % Channels
    if isfield(in,'channel')
        if ~isempty(in.channel)
            Chn = {ChanMat.Channel.Name};
            has.channel(iInput) = ruleFcn(ismember(upper(in.channel),upper(Chn)));
        end
    end
    
    % Events
    if isfield(in,'event')
        if ~isempty(in.event)
            Evt = {DataMat.F.events.label};
            has.event(iInput) = ruleFcn(ismember(upper(in.event),upper(Evt)));
        end
    end

end

end

% -------------------------------------------------------------------------
function Report(sInputs, sProcess, FeatName, in, has)
% Subroutine of Process.

SubjectNames = {sInputs.SubjectName};

switch sProcess.options.rule.Value
    case 1 % Files have (or not) ALL
        strRule = 'ALL';
    case 2 % Files have (or not) ANY
        strRule = 'ANY';
end


if ~isempty(in)
    strMsg = [];
    strMsg = [strMsg, sprintf('Following %s(s) were specified:',FeatName)];
    strMsg = [strMsg, sprintf('\n')];
    strMsg = [strMsg, sprintf('\n - %s',in{:})];
    strMsg = [strMsg, sprintf('\n\n')];
    strMsg = [strMsg, sprintf('%d file(s) HAVE [%s] specified %s(s):',numel(find(has)),strRule,FeatName)];
    strMsg = [strMsg, sprintf('\n')];
    strMsg = [strMsg, sprintf('\n - %s',SubjectNames{has})];
    strMsg = [strMsg, sprintf('\n\n')];
    strMsg = [strMsg, sprintf('%d file(s) HAVE NOT [%s] specified %s(s):',numel(find(~has)),strRule,FeatName)];
    strMsg = [strMsg, sprintf('\n')];
    strMsg = [strMsg, sprintf('\n - %s',SubjectNames{~has})];
    % Using Warning for Report Viewer to open after process
    bst_report('Warning', sProcess, sInputs, strMsg);
end

end
