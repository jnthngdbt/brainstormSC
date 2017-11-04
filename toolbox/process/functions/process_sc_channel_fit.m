function varargout = process_sc_channel_fit( varargin )
% process_sc_channel_fit: Fit a custom channel montage with a template
% channel file from Brainstorm. Fit is based on channel labels, so order of
% labels will adapt to file's montage.
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % ===== DEFAULTS CHANNEL FILE TEMPLATES =====
    strList = {''};
    % Get registered Brainstorm EEG defaults
    bstDefaults = bst_get('EegDefaults');
    % Build a list of strings representing all the defaults
    for iGroup = 1:length(bstDefaults)
        for iDef = 1:length(bstDefaults(iGroup).contents)
            strList{end+1} = [bstDefaults(iGroup).name ': ' bstDefaults(iGroup).contents(iDef).name];
        end
    end
    
    % ===== PROCESS =====
    % Description the process
    sProcess.Comment     = 'CHANNEL > Fit channel file based on labels';
    sProcess.FileTag     = '';
    sProcess.Description = 'Fit a custom channel montage from a template channel file from Brainstorm. Fit based on channel labels.';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(345);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    % Option: Default channel files
    sProcess.options.template.Comment = 'Template channel file to fit:';
    sProcess.options.template.Type    = 'combobox';
    sProcess.options.template.Value   = {1, strList};
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
        
% ===== GET CHANNEL TEMPLATE TO FIT =====
% Get registered Brainstorm EEG defaults
bstDefaults = bst_get('EegDefaults');
% Get default channel file
iSel   = sProcess.options.template.Value{1};
strDef = sProcess.options.template.Value{2}{iSel};
% If there is something selected
if ~isempty(strDef)
    % Format: "group: name"
    cDef   = strtrim(str_split(strDef, ':'));
    % Find the selected group in the list 
    iGroup = find(strcmpi(cDef{1}, {bstDefaults.name}));
    % If group is found
    if ~isempty(iGroup)
        % Find the selected default in the list 
        iDef = find(strcmpi(cDef{2}, {bstDefaults(iGroup).contents.name}));
        % If default was found
        if ~isempty(iDef)
            sDef = bstDefaults(iGroup).contents(iDef);
        end
    end
end

% === COMPUTE FIT ON ALL FILES ===
bst_progress('start', 'Channel fit based on labels', 'Processing files...', 0, numel(sInputs));
for iInput=1:numel(sInputs)
    bst_progress('set', iInput);
    % Get file's channel file
    sChan = in_bst_channel(sInputs(iInput).ChannelFile);
    
    % Compute the fit
    sChan = Fit(sChan, sDef);  

    % ===== SET CHANNEL FILE =====
    % Get channel studies
    [tmp, iChanStudies] = bst_get('ChannelForStudy', [sInputs(iInput).iStudy]);
    iChanStudies = unique(iChanStudies);
    % Set options
    ChannelAlign = 0;% 2 * double(sProcess.options.channelalign.Value);
    ChannelReplace = 2;
    UpdateData = 1;
    % Update channel file in database
    db_set_channel(iChanStudies, sChan, ChannelReplace, ChannelAlign, UpdateData);
end

% Return all the files in input
OutputFiles = {sInputs.FileName};

end

%% ===== COMPUTE =====
function sChan = Fit(sChan,sDef)
% Fit a custom channel montage with a template channel file from 
% Brainstorm, mainly to obtain sensors location when unavailbable. Fit is 
% based on channel labels. Order of input channels doesn't change.
% 
%   sChan = Fit(sChan,sDef)
% 
% Input 'sChan' is the channel file structure as obtained from 
% <in_bst_channel>. 'sDef' is one of the structures of default channel 
% files as obtained from <bst_get('EegDefaults')>, corresponding to the 
% channel file to fit. The struct has 2 fields: 'fullpath' and 'name'

% Load template channel file to fit
sChanTemp = load(sDef.fullpath);

% Custom channel file labels
labels = {sChan.Channel.Name};
% Template channel file labels
tempLabels = {sChanTemp.Channel.Name};

% Scan through each label, find a fit in template
for ii=1:numel(labels)
    % Fitting labels
    iFit = strcmpi(labels{ii},tempLabels);
    if any(iFit) % If fit found, copy specific channel struct
        sChan.Channel(ii) = orderfields(sChanTemp.Channel(iFit),sChan.Channel(ii));
    else
        if sum(sChan.Channel(ii).Loc)==0
            sChan.Channel(ii).Loc = 0.001*randn(3,1);
        end
        if      ~isempty(strfind(labels{ii},'ECG')), sChan.Channel(ii).Type = 'ECG';
        elseif  ~isempty(strfind(labels{ii},'EKG')), sChan.Channel(ii).Type = 'ECG';
        elseif  ~isempty(strfind(labels{ii},'EOG')), sChan.Channel(ii).Type = 'EOG';
        elseif  ~isempty(strfind(labels{ii},'EMG')), sChan.Channel(ii).Type = 'EMG';
        else    sChan.Channel(ii).Type = 'MISC';
        end
    end
end

end