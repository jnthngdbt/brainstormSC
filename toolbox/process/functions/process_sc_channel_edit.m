function varargout = process_sc_channel_edit( varargin )
% process_sc_channel_edit: Edit informations in channel file in a script
% way
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    
    % ===== PROCESS =====
    % Description the process
    sProcess.Comment     = 'CHANNEL > Edit channel file';
    sProcess.FileTag     = '';
    sProcess.Description = 'Edit informations in channel file in a script way';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(335);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    %
    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all)';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'EEG';
    %
    sProcess.options.prefix.Comment = 'Add prefix';
    sProcess.options.prefix.Type    = 'text';
    sProcess.options.prefix.Value   = [];
    %
    sProcess.options.suffix.Comment = 'Add suffix';
    sProcess.options.suffix.Type    = 'text';
    sProcess.options.suffix.Value   = [];
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
  
% Get input options
prefix = sProcess.options.prefix.Value;
suffix = sProcess.options.suffix.Value;
sensor = sProcess.options.sensortypes.Value;

for iInput=1:numel(sInputs)
    % Get channel struct
    sChan = in_bst_channel(sInputs(iInput).ChannelFile);
    % Get selected channel indexes
    [iChannels, Comment] = channel_find(sChan.Channel, sensor);
    
    % Add prefix
    if ~isempty(prefix)
        for iChan = iChannels
            sChan.Channel(iChan).Name = [prefix,sChan.Channel(iChan).Name];
        end
    end
    
    % Add suffix
    if ~isempty(suffix)
        for iChan = iChannels
            sChan.Channel(iChan).Name = [sChan.Channel(iChan).Name,suffix];
        end
    end
    
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

