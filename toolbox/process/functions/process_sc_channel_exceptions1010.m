function varargout = process_sc_channel_exceptions1010( varargin )
% process_sc_channel_exceptions1010: Rename EEG 10-10 and 10-20 label systems to
% deal with exceptions:
% 
%         | 10-20 | 10-10 |
% --------+-------+-------+
% 'T7/T3' | 'T3'  | 'T7'  |
% 'T8/T4' | 'T4'  | 'T8'  |
% 'P7/T5' | 'T5'  | 'P7'  |
% 'P8/T6' | 'T6'  | 'P8'  |
% 
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>    
    % ===== PROCESS =====
    % Description the process
    sProcess.Comment     = 'CHANNEL > Rename label exceptions in 10-10/10-20 EEG systems';
    sProcess.FileTag     = '';
    sProcess.Description = 'Rename label exceptions in 10-10/10-20 EEG systems';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(340);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
%     % Separator
%     sProcess.options.separator.Type = 'separator';
%     sProcess.options.separator.Comment = ' ';
    %
    sProcess.options.system.Comment = {...
        'Use both systems labels (T7/T3, T8/T4, P7/T5, P8/T6)', ...
        'Use 10-20 system labels (T3, T4, T5, T6)', ...
        'Use 10-10 system labels (T7, T8, P7, P8)'};
    sProcess.options.system.Type    = 'radio';
    sProcess.options.system.Value   = 1;
    % === 
    sProcess.options.label2.Comment = [ ...
        '<HTML><BR>- Channel files "10-10 65" and "10-20 19" use both systems' ...
        '<BR>  - Channel files "BioSemi *" and "BrainProducts *" use 10-10 system', ...
        ];
    sProcess.options.label2.Type    = 'label';

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
        
% Wich system to use
switch sProcess.options.system.Value
    case 1
        labelSystem = 'both';
    case 2
        labelSystem = '10-20';
    case 3
        labelSystem = '10-10';
end

% === COMPUTE ON ALL FILES ===
for iInput=1:numel(sInputs)
    
    % Get file's channel file struct
    sChan = in_bst_channel(sInputs(iInput).ChannelFile);
    % Channels labels
    labels = {sChan.Channel.Name};
     
    % Process labels
    labels = RenameExceptionLabels(labels,labelSystem);     
    
    % Reput labels in channel struct
    for iChan=1:numel(labels)
        sChan.Channel(iChan).Name = labels{iChan};
    end
    
    % ===== UPDATE CHANNEL FILE =====
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

% ==
function labels = RenameExceptionLabels(labels,labelSystem)

switch lower(labelSystem)
    case 'both'
        labels(strcmpi('T8',labels)) = {'T8/T4'};
        labels(strcmpi('T4',labels)) = {'T8/T4'};
        labels(strcmpi('T7',labels)) = {'T7/T3'};
        labels(strcmpi('T3',labels)) = {'T7/T3'};
        labels(strcmpi('P8',labels)) = {'P8/T6'};
        labels(strcmpi('T6',labels)) = {'P8/T6'};
        labels(strcmpi('P7',labels)) = {'P7/T5'};
        labels(strcmpi('T5',labels)) = {'P7/T5'};
    case '10-20'
        labels(strcmpi('T8',labels))    = {'T4'};
        labels(strcmpi('T8/T4',labels)) = {'T4'};
        labels(strcmpi('T7',labels))    = {'T3'};
        labels(strcmpi('T7/T3',labels)) = {'T3'};
        labels(strcmpi('P8',labels))    = {'T6'};
        labels(strcmpi('P8/T6',labels)) = {'T6'};
        labels(strcmpi('P7',labels))    = {'T5'};
        labels(strcmpi('P7/T5',labels)) = {'T5'};
    case '10-10'
        labels(strcmpi('T4',labels))    = {'T8'};
        labels(strcmpi('T8/T4',labels)) = {'T8'};
        labels(strcmpi('T3',labels))    = {'T7'};
        labels(strcmpi('T7/T3',labels)) = {'T7'};
        labels(strcmpi('T6',labels))    = {'P8'};
        labels(strcmpi('P8/T6',labels)) = {'P8'};
        labels(strcmpi('T5',labels))    = {'P7'};
        labels(strcmpi('P7/T5',labels)) = {'P7'};
end

end


