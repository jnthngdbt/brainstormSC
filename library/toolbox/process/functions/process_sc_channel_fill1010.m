function varargout = process_sc_channel_fill1010( varargin )
% process_sc_channel_fill1010: Rename EGI GSN257 hat electrodes
% corresponding to 10-20 EEG system with interpretable 10-20 labels (ex:
% E37->Fp1).
% 
% External calls:
% 
%       labels = process_sc_channel_fill1010('FillGSN257With1010',labels)
% 
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>    
    % ===== PROCESS =====
    % Description the process
    sProcess.Comment     = 'CHANNEL > HYDROCEL GSN257 > Rename some labels to EEG 10-20';
    sProcess.FileTag     = '';
    sProcess.Description = 'Rename some labels to EEG 10-20 (ex: E37->Fp1)';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(355);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    % === 
    sProcess.options.label1.Comment = '<HTML>Rename some labels to EEG 10-20 (ex: E37->Fp1)';
    sProcess.options.label1.Type    = 'label';
    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
    % 
    sProcess.options.simulateaux.Comment = 'Rename some EEG to simulate EOG and EMG';
    sProcess.options.simulateaux.Type    = 'checkbox';
    sProcess.options.simulateaux.Value   = 0;
    % 
    sProcess.options.roc.Comment = 'Use following for ROC (right occular): ';
    sProcess.options.roc.Type    = 'text';
    sProcess.options.roc.Value   = 'E1';
    % 
    sProcess.options.loc.Comment = 'Use following for LOC (left occular): ';
    sProcess.options.loc.Type    = 'text';
    sProcess.options.loc.Value   = 'E241';
    % 
    sProcess.options.emgr.Comment = 'Use following for EMGR (right muscular): ';
    sProcess.options.emgr.Type    = 'text';
    sProcess.options.emgr.Value   = 'E240';
    % 
    sProcess.options.emgl.Comment = 'Use following for EMGL (left muscular): ';
    sProcess.options.emgl.Type    = 'text';
    sProcess.options.emgl.Value   = 'E243';
    % === 
    sProcess.options.label2.Comment = '<HTML>Suggested: E1->ROC, E241->LOC, E240->EMGR, E243->EMGL';
    sProcess.options.label2.Type    = 'label';

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
        
% === COMPUTE ON ALL FILES ===
for iInput=1:numel(sInputs)
    
    % Get file's channel file struct
    sChan = in_bst_channel(sInputs(iInput).ChannelFile);
    % Channels labels
    labels = {sChan.Channel.Name};
     
    % Process labels
    labels = Normalize1010Labels(labels); 
    labels = FillGSN257With1010(labels);    % EGI GSN257 Hat
    
    % Simulate auxilaries (EOG, EMG)
    if sProcess.options.simulateaux.Value>0
        labels(strcmpi(sProcess.options.roc.Value,labels)) = {'ROC'};
        labels(strcmpi(sProcess.options.loc.Value,labels)) = {'LOC'};
        labels(strcmpi(sProcess.options.emgl.Value,labels)) = {'EMGL'};
        labels(strcmpi(sProcess.options.emgr.Value,labels)) = {'EMGR'};
    end
    
    % Reput labels in channel struct
    for iChan=1:numel(labels)
        sChan.Channel(iChan).Name = labels{iChan};
    end
    
    % Make sure Cz (VREF) is of type EEG (not fiducial)
    iCz = strcmpi({sChan.Channel.Name},'Cz');
    sChan.Channel(iCz).Type = 'EEG';
    
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
function labels = Normalize1010Labels(labels)
% 10-10
labels(strcmpi('T8',labels)) = {'T8/T4'};
labels(strcmpi('T4',labels)) = {'T8/T4'};
labels(strcmpi('T7',labels)) = {'T7/T3'};
labels(strcmpi('T3',labels)) = {'T7/T3'};
labels(strcmpi('P8',labels)) = {'P8/T6'};
labels(strcmpi('T6',labels)) = {'P8/T6'};
labels(strcmpi('P7',labels)) = {'P7/T5'};
labels(strcmpi('T5',labels)) = {'P7/T5'};

end

% ==
function labels = FillGSN257With1010(labels)

% GSN 257
labels(strcmpi('E37',labels)) = {'Fp1'};
labels(strcmpi('E18',labels)) = {'Fp2'};
labels(strcmpi('E36',labels)) = {'F3'};
labels(strcmpi('E21',labels)) = {'Fz'};
labels(strcmpi('E224',labels)) = {'F4'};
labels(strcmpi('E59',labels)) = {'C3'};
labels(strcmpi('E257',labels)) = {'Cz'};
labels(strcmpi('E183',labels)) = {'C4'};
labels(strcmpi('E87',labels)) = {'P3'};
labels(strcmpi('E101',labels)) = {'Pz'};
labels(strcmpi('E153',labels)) = {'P4'};
labels(strcmpi('E116',labels)) = {'O1'};
labels(strcmpi('E126',labels)) = {'Oz'};
labels(strcmpi('E150',labels)) = {'O2'};
labels(strcmpi('E2',labels)) = {'F8'};
labels(strcmpi('E47',labels)) = {'F7'};
labels(strcmpi('E94',labels)) = {'TP9'};
labels(strcmpi('E190',labels)) = {'TP10'};
labels(strcmpi('E202',labels)) = {'T8/T4'};
labels(strcmpi('E69',labels)) = {'T7/T3'};
labels(strcmpi('E170',labels)) = {'P8/T6'};
labels(strcmpi('E96',labels)) = {'P7/T5'};

end

