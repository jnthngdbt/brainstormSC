function varargout = process_sc_reject_trials_evt( varargin )
% process_sc_reject_trials_evt: Reject trials based on events

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'ARTIFACTS > Reject trials based on events';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(545);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    % 
    sProcess.options.event.Comment = 'Event name(s):';
    sProcess.options.event.Type    = 'text';
    sProcess.options.event.Value   = '';
    % 
    sProcess.options.coverage.Comment = 'Minimum of data affected (only for extended artifact events):';
    sProcess.options.coverage.Type    = 'value';
    sProcess.options.coverage.Value   = {10,'%',0};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

% === INPUT PARAMETERS ===

% Get event names specified
artName = [];
if ~isempty(sProcess.options.event.Value)
    artName = textscan(sProcess.options.event.Value,'%s','delimiter',',');
    artName = artName{1};
end

% === COMPUTE ===

OutputFiles = [];
Ni = numel(sInputs);
isBadTrials = false(1,Ni);
bst_progress('start','Process','Looking for artifact in files',0,Ni);
for iInput=1:Ni
    bst_progress('set',iInput);
    
    % Load data structs
    sData = in_bst_data(sInputs(iInput).FileName);
    % Check if data has events
    if ~isfield(sData,'Events'), continue; end
    if isempty(sData.Events),    continue; end
    
    % Get artifact events in data
    iArt = find(ismember({sData.Events.label},artName));
    isArt = false(1,numel(sData.Time));
    if isempty(iArt),    continue; end
    for i=iArt
        % Time position of artifacts in data
        times = [sData.Events(iArt).times];
        [Np,Ne] = size(times);
        isExtended = Np==2;
        if isExtended
            % Flag artifacted data time points
            for k=1:Ne
                isArt = isArt | (sData.Time >= times(1,k) & sData.Time <= times(2,k));
            end
        end 
    end
    % 
    if isExtended 
        % Mark trial as bad if artifacts covers more than a percentage
        isBadTrials(iInput) = numel(find(isArt))/numel(isArt) >= sProcess.options.coverage.Value{1}/100;
    else
        % Mark trial as bad if any artifact event is in the data
        isBadTrials(iInput) = ~isempty(iArt);
    end

end

% Record bad trials in study (extracted from PROCESS_DETECTBAD)
iBadTrials = find(isBadTrials);
if ~isempty(iBadTrials)
    isBad = true;
    process_detectbad('SetTrialStatus',{sInputs(iBadTrials).FileName}, isBad);
    process_detectbad('SetTrialStatus',{sInputs(~iBadTrials).FileName}, ~isBad);
    bst_report('Info', sProcess, sInputs, sprintf('Epochs tested: %d - Bad epochs: %d (%d%%)', length(sInputs), length(iBadTrials), round(nnz(iBadTrials)/length(sInputs)*100)));
end
% Return only good trials
iGoodTrials = setdiff(1:length(sInputs), iBadTrials);
OutputFiles = {sInputs(iGoodTrials).FileName};

end


% % First pass to check if all files have what is wanted
% notok = 0;
% for iInput=1:numel(sInputs)
%     % Load the raw file descriptor
%     DataMat = in_bst_data(sInputs(iInput).FileName);
%     sFile = DataMat.F;
%        
%     eventi = {sFile.events.label};
%     idx = find(~ismember(lower(eventName),lower(eventi)));
%     if ~isempty(idx), notok = 1; end
% 
% end
% 
% if notok, error('All files do not have all that is wanted'); end



