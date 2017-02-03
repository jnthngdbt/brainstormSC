function varargout = process_sc_report_venn( varargin )
% process_sc_report_venn: Reports what is common and different
% (channels, events, ...) in a group of datasets (RAW files) using Venn
% analysis
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'REPORT > Database Venn analysis';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(245);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % 
    sProcess.options.label1.Comment = '<HTML>Reports channels and events that are common (intersection)<BR> and not (union-intersection) across files.';
    sProcess.options.label1.Type    = 'label';

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

bst_progress('start', 'Processing', 'Extracting "Venn" information from files...', 0, numel(sInputs));
for iInput=1:numel(sInputs)
    bst_progress('set', iInput);
    
    % Get data and channel structs
    DataMat = in_bst_data(sInputs(iInput).FileName,sInputs(iInput).ChannelFile);
    ChanMat = in_bst_channel(sInputs(iInput).ChannelFile);
            
    % Accumulate Venn sets (union and intersect)
    if iInput==1 % First pass
        % Channels
        unChn = {ChanMat.Channel.Name};
        intChn = unChn;
        
        % Events
        unEvt = {DataMat.F.events.label};
        intEvt = unEvt;
        
    else
        % Channels
        Chn = {ChanMat.Channel.Name};
        [unChn, intChn] = Venn(Chn, unChn, intChn);
        
        % Events
        Evt = {DataMat.F.events.label};
        [unEvt, intEvt] = Venn(Evt, unEvt, intEvt);
        
    end
    
end

% Set differences
diffChn = setdiff(unChn,intChn);
diffEvt = setdiff(unEvt,intEvt);

% % % % % % Determine the number of files than differs
% Interesting, but may be take too long to span all setdifs, since have to
% scan all input files every time
% % % % % ndiffChn = zeros(1,diffChn);
% % % % % for i=1:numel(diffChn)
% % % % %     [has] = Has(sInputs, in);
% % % % % end

% === REPORT ===

Report(sProcess, sInputs, 'channel', intChn, diffChn);
Report(sProcess, sInputs, 'event', intEvt, diffEvt);

% bst_report('Open','last')
% bst_report('Open','next')
% bst_report('Open','loaded')

% === EXIT ===
OutputFiles = {sInputs.FileName}; % Return input file names

end

% -------------------------------------------------------------------------
function [unVal, intVal] = Venn(Val, unVal, intVal)

if isempty(Val)
    intVal = [];
    return;
end
% Union
[c,iA,iB] = union(upper(unVal),upper(Val));
unVal = [unVal(iA),Val(iB)]; % To keep original case
% Intersect
[c,iA,iB] = intersect(upper(intVal),upper(Val));
intVal = intVal(iA); % To keep original case

end

% -------------------------------------------------------------------------
function Report(sProcess, sInputs, FeatName, Int, Diff)

strMsg = [];
strMsg = [strMsg, sprintf('Following %s(s) are common to all files (intersection):',FeatName)];
strMsg = [strMsg, sprintf('\n')];
strMsg = [strMsg, sprintf('\n - %s',Int{:})];
% bst_report('Info', sProcess, sInputs, strMsg);

strMsg = [strMsg, sprintf('\n\n')];
strMsg = [strMsg, sprintf('Following %s(s) are NOT common to all files (union/intersection):',FeatName)];
strMsg = [strMsg, sprintf('\n')];
strMsg = [strMsg, sprintf('\n - %s',Diff{:})];

% Using Warning for Report Viewer to open after process
bst_report('Warning', sProcess, sInputs, strMsg);  

end
