function varargout = process_sc_import_data_events( varargin )
% PROCESS_SC_IMPORT_DATA_EVENTS: Import data from file using event markers
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'IMPORT > DATA > From events';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(75);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw',  'data'};
    sProcess.OutputTypes = {'data', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    % 
    sProcess.options.event.Comment = 'Event name(s):';
    sProcess.options.event.Type    = 'text';
    sProcess.options.event.Value   = '';
    % 
    sProcess.options.condition.Comment = 'Group into one condition having this name (Empty: one condition per event type):';
    sProcess.options.condition.Type    = 'text';
    sProcess.options.condition.Value   = '';
    % Separator -----------------------------------
    sProcess.options.sep1.Type = 'separator';
    sProcess.options.sep1.Comment = ' ';
    %
    sProcess.options.useext.Comment = 'Use events extension as epoch time (otherwise, specify below)';
    sProcess.options.useext.Type    = 'checkbox';
    sProcess.options.useext.Value   = 0;
    % Epoch time
    sProcess.options.epochtime.Comment = 'Epoch time: ';
    sProcess.options.epochtime.Type    = 'range';
    sProcess.options.epochtime.Value   = {[0, 2], 's', []};
    % 
    sProcess.options.alignment.Comment = {'Time origin is first sample', 'Time origin is middle sample','Time origin is last sample'};
    sProcess.options.alignment.Type    = 'radio';
    sProcess.options.alignment.Value   = 1;
    % Separator -----------------------------------
    sProcess.options.sep2.Type = 'separator';
    sProcess.options.sep2.Comment = ' ';
    % === 
    sProcess.options.label1.Comment = '<HTML>Brainstorm import function options: ';
    sProcess.options.label1.Type    = 'label';
    % Use SSP
    sProcess.options.usessp.Comment = 'Use SSP projectors (e.g. for artifacts rejection)';
    sProcess.options.usessp.Type    = 'checkbox';
    sProcess.options.usessp.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

OutputFiles = {};

% === INPUT PARAMETERS ===

% Event names
EvtNames = strtrim(str_split(sProcess.options.event.Value, ',;'));
if isempty(EvtNames)
    bst_report('Error', sProcess, [], 'No events to import.');
    return;
end

% Use events extension?
useExtension = sProcess.options.useext.Value;

% Time window to import
epochtime = sProcess.options.epochtime.Value{1};

% Import options
ImportOptions = db_template('ImportOptions');
ImportOptions.ImportMode        = 'Event';
ImportOptions.TimeRange         = [];%sProcess.options.timewindow.Value{1};
ImportOptions.UseEvents         = 1;
ImportOptions.EventsTimeRange   = epochtime;
ImportOptions.iEpochs           = 1;
ImportOptions.SplitRaw          = 0;
ImportOptions.UseCtfComp        = 0; %sProcess.options.usectfcomp.Value;
ImportOptions.UseSsp            = sProcess.options.usessp.Value;
% % % % % % % % ImportOptions.CreateConditions  = 0; % We create custom condition, not based on event type
ImportOptions.ChannelReplace    = 0;
ImportOptions.ChannelAlign      = 0;
ImportOptions.IgnoreShortEpochs = 0;
ImportOptions.EventsMode        = 'ignore';
ImportOptions.EventsTrackMode   = 'value';
ImportOptions.DisplayMessages   = 0;

% Condition
if isempty(sProcess.options.condition.Value)
    Condition = [];
    ImportOptions.CreateConditions  = 1;
else
    Condition = sProcess.options.condition.Value;
    ImportOptions.CreateConditions  = 0;
end

% === COMPUTE ===

Ni = numel(sInputs);

% Initialize progress bar
if bst_progress('isVisible'), bst_progress('set', bst_progress('get')); end
bst_progress('start', 'Processing', 'Importing from events...', 0, Ni);

% Import files
FileNames = {sInputs.FileName};
for iFile = 1:Ni
    bst_progress('set', iFile);
% % % % %     % Condition name
% % % % %     if isempty(sProcess.options.condition.Value)
% % % % %         Condition = sInputs(iFile).Condition;
% % % % %         if isempty(Condition)
% % % % %             Condition = EvtNames{1};
% % % % %         end
% % % % %     else
% % % % %         Condition = sProcess.options.condition.Value;
% % % % %     end
    
    isRaw = strcmpi(sInputs(iFile).FileType, 'raw');
    if isRaw
        FileMat    = in_bst_data(sInputs(iFile).FileName, [], 'F');
        sFile      = FileMat.F;
        FileFormat = sFile.format;
    else
        FileFormat = 'BST-DATA';
        sFile = in_fopen(sInputs(iFile).FileName, FileFormat);
    end
    
%     FileMat = in_bst_data(FileNames{iFile}, [], 'F');
%     sFile   = FileMat.F;

    % ===== IMPORT FILES =====
    SubjectName = sInputs(iFile).SubjectName;
    % Get subject
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    % Create subject if it does not exist yet
    if isempty(sSubject)
        [sSubject, iSubject] = db_add_subject(SubjectName);
    end
    % Define output study
    if ~isempty(Condition)
        % Get condition asked by user
        [sStudy, iStudy] = bst_get('StudyWithCondition', fullfile(SubjectName, Condition));
        % Condition does not exist: create it
        if isempty(sStudy)
            iStudy = db_add_condition(SubjectName, Condition, 1);
            if isempty(iStudy)
                bst_report('Error', sProcess, sInputs, ['Cannot create condition : "' fullfile(SubjectName, Condition) '"']);
                return;
            end
        end
    else
        iStudy = [];
    end


    % No events in file
    if isempty(sFile.events)
        bst_report('Error', sProcess, [], ['No events in file: ' 10 FileNames{iFile}]);
        continue;
    end
    % Initialize events structure
    events = repmat(sFile.events, 0);
    % Get selected events
    for iSelEvt = 1:length(EvtNames)
        % Find input event in file
        iEvt = find(strcmpi(EvtNames{iSelEvt}, {sFile.events.label}));
        if isempty(iEvt)
            bst_report('Warning', sProcess, [], ['Event "' EvtNames{iSelEvt} '" does not exist in file: ' 10 FileNames{iFile}]);
            continue;
        end
        newEvt = sFile.events(iEvt);
        % Add to the list of events to import
        events(end+1) = newEvt;
    end
    % No events to import found in the file
    if isempty(events)
        bst_report('Error', sProcess, [], ['No events to import found in file: ' 10 FileNames{iFile}]);
        continue;
    end
    % Define epoch time to import
    if ~useExtension
        for iEvt=1:numel(events)
            switch sProcess.options.alignment.Value
                case 1
                    events(iEvt).samples = events(iEvt).samples(1,:);
                    events(iEvt).times   = events(iEvt).times(1,:);
                case 2
                    if size(events.samples,1)==2
                        events(iEvt).samples = round(mean(events(iEvt).samples));
                        events(iEvt).times   = mean(events(iEvt).times);
                    else
                        events(iEvt).samples = events(iEvt).samples(1,:);
                        events(iEvt).times   = events(iEvt).times(1,:);
                    end
                case 3
                    if size(events.samples,1)==2
                        events(iEvt).samples = events(iEvt).samples(2,:);
                        events(iEvt).times   = events(iEvt).times(2,:);
                    else
                        events(iEvt).samples = events(iEvt).samples(1,:);
                        events(iEvt).times   = events(iEvt).times(1,:);
                    end
            end
        end
    else
        for iEvt=1:numel(events)
            if size(events(iEvt).samples,1)==1
                bst_report('Error', sProcess, [], 'It was specified to use events extension (duration), but events are not extended.');
                continue;
            end
        end
    end
    ImportOptions.events = events;
    % Import file
    OutputFiles = cat(2, OutputFiles, import_data(sFile, FileFormat, iStudy, iSubject, ImportOptions));
end
% Report number of files generated
if ~isempty(OutputFiles)
    bst_report('Info', sProcess, sInputs, sprintf('%d epochs imported.', length(OutputFiles)));
else
    bst_report('Error', sProcess, sInputs, 'Nothing imported from those files.');
end

end

% for iInput=1:Ni
%     
%     % Load the raw file descriptor
%     DataMat = in_bst_data(sInputs(iInput).FileName);
%     sFile = DataMat.F;
%     % Sampling frequency
%     Fs = 1/(DataMat.Time(2)-DataMat.Time(1));
%     % File events
%     FileEvt = {sFile.events.label};
%     
%     Samples = [];
%     for k=1:numel(eventName)
%         % Get event
%         iEvt = find(ismember(lower(FileEvt),lower(eventName{k})));
%         if isempty(iEvt)
%             bst_report('Warning', sProcess, sInputs(iInput), ['Event [' eventName{k} '] unavailable']);
%             continue;
%         end
%         sEvt = sFile.events(iEvt);
% 
%         Samples = [Samples, sEvt.samples(1,:)];
%     end
%     
%     % If no events were found, continue to next file
%     if isempty(Samples)
%         bst_report('Warning', sProcess, sInputs(iInput), 'No data imported');
%         continue;
%     end
%     
%     % Get samples bounds
%     Ne = size(Samples,2);
%     Samples = sort(Samples,2,'ascend');
%     SamplesBounds = zeros(2,Ne);
%     SamplesBounds(1,:) = Samples(1,:) + round(epochtime(1)*Fs);
%     SamplesBounds(2,:) = Samples(1,:) + round(epochtime(2)*Fs);
%     
%     % Extract data matrices for each event
%     bst_progress('start','Import process',sprintf('Import %d events of file %d/%d',Ne,iInput,Ni),0,Ne);
%     for i=1:Ne
%         bst_progress('set',i);
%         
% %         [F, TimeVector] = in_fread(sFile, 1, SamplesBounds, iChannels);
%         [F, TimeVector] = in_fread(sFile, 1, SamplesBounds(:,i));
%         [Nc,Nt] = size(F);
%         TimeVector = linspace(epochtime(1),epochtime(2),Nt);
%         
%         DataMat = db_template('datamat');
%         DataMat.F = F;
%         DataMat.Comment = sprintf('%s (#%d)',condName,i);
%         DataMat.ChannelFlag = ones(Nc,1);
%         DataMat.Time = TimeVector;
%         DataMat.DataType = 'data';
%         DataMat.Device = sFile.device;
%         DataMat.nAvg = 1;
% %         DataMat.Events = struct([]);
% 
%         % --------
%         % Add events included in time window
%         
%         % Initialize events structure
%         events = repmat(sFile.events, 0);
%         if ~isempty(sFile.events)
%             % 
%             for iEvt = 1:length(FileEvt)
%                 [iInput,i,iEvt]
%                 newEvt = sFile.events(iEvt);
%                 isExtended = size(newEvt.samples,1)==2;
%                 % Find events that are in time window
%                 iOcc = newEvt.samples(1,:) < SamplesBounds(2,i);
%                 if isExtended
%                     iOcc = iOcc & (newEvt.samples(2,:) > SamplesBounds(1,i));
%                 else
%                     iOcc = iOcc & (newEvt.samples(1,:) > SamplesBounds(1,i));
%                 end
%                 iOcc = find(iOcc);
%                 % No occurrence for this event: skip to next event
%                 if isempty(iOcc)
%                     continue;
%                 end
%                 % Get the selected occurrences
%                 newEvt.epochs  = newEvt.epochs(iOcc);
%                 newEvt.samples = newEvt.samples(:, iOcc);
%                 newEvt.times   = newEvt.times(:, iOcc);
%                 % Add to the list of events to import
%                 events(end+1) = newEvt;
%             end
%         end
%         DataMat.Events = events;
%         % --------
%         
%         SubjectName = sFile.comment;
%         Condition = condName;
%         % Get subject
%         [sSubject, iSubject] = bst_get('Subject', SubjectName);
%         % Create subject if it does not exist yet
%         if isempty(sSubject)
%             [sSubject, iSubject] = db_add_subject(SubjectName);
%         end
%         % Define output study
%         % Get condition asked by user
%         [sStudy, iStudy] = bst_get('StudyWithCondition', fullfile(SubjectName, Condition));
%         % Condition does not exist: create it
%         if isempty(sStudy)
%             iStudy = db_add_condition(SubjectName, Condition, 1);
%             if isempty(iStudy)
%                 bst_report('Error', sProcess, sInputs, ['Cannot create condition : "' fullfile(SubjectName, Condition) '"']);
%                 return;
%             end
%         end
%         sStudy = bst_get('Study',iStudy);
% 
%         % Output filename: add file tag
%         sProtocol = bst_get('ProtocolInfo');
%         OutputFile = fullfile(sProtocol.STUDIES, SubjectName, sStudy.Name, sprintf('data_%s_trial.mat',sStudy.Name));
%         OutputFile = file_unique(OutputFile);
%         % Save file
%         bst_save(OutputFile, DataMat, 'v6');
%         % Add file to database structure
%         db_add_data(iStudy, OutputFile, DataMat);
%         OutputFiles = [OutputFiles; {OutputFile}];
%     end
% 
% end




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



