function varargout = process_sc_concat( varargin )
% PROCESS_CONCAT: Concatenate several files using the time information from
% the first file. The difference with built-in BST concat process is that
% it handles timefreq data (July 2013). But most code has been copy-pasted
% from BST function!
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'FILE > Concatenate (time dimension)';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(845);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data','timefreq'};
    sProcess.OutputTypes = {'data','timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 2;
    sProcess.isSeparator = 1;
    % Definition of the options
    % === AVERAGE TYPE
    sProcess.options.avgtype.Comment = {'Everything', 'By subject', 'By condition (subject-level)', 'By condition (all)', 'By trial group (subject-level)', 'By trial group (all)'};
    sProcess.options.avgtype.Type    = 'radio';
    sProcess.options.avgtype.Value   = 1;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    % Initialize returned list
    OutputFiles = {};
    % Get current progressbar position
    if bst_progress('isVisible')
        curProgress = bst_progress('get');
    else
        curProgress = [];
    end
    % Group files in different ways: by subject, by condition, or all together
    switch(sProcess.options.avgtype.Value)
        % === EVERYTHING ===
        case 1  
            OutputFiles = Concatenate(sProcess, sInputs);
            
        % === BY SUBJECT ===
        case 2
            % Process each subject independently
            uniqueSubj = unique({sInputs.SubjectFile});
            for i = 1:length(uniqueSubj)
                % Set progress bar at the same level for each loop
                if ~isempty(curProgress)
                    bst_progress('set', curProgress);
                end
                % Get all the files for condition #i
                iInputSubj = find(strcmpi(uniqueSubj{i}, {sInputs.SubjectFile}));
                % Process the average of condition #i
                tmpOutFiles = Concatenate(sProcess, sInputs(iInputSubj));
                % Save the results
                OutputFiles = cat(2, OutputFiles, tmpOutFiles);
            end
            
        % === BY CONDITION ===
        case {3,4}
            % Subject average
            if (sProcess.options.avgtype.Value == 3)
                inputCondPath = {};
                for iInput = 1:length(sInputs)
                    inputCondPath{iInput} = [sInputs(iInput).SubjectName, '/', sInputs(iInput).Condition];
                end
            % Grand average
            else
                inputCondPath = {sInputs.Condition};
            end
            % Process each condition independently
            uniqueCond = unique(inputCondPath);
            for i = 1:length(uniqueCond)
                % Set progress bar at the same level for each loop
                if ~isempty(curProgress)
                    bst_progress('set', curProgress);
                end
                % Get all the files for condition #i
                iInputCond = find(strcmpi(uniqueCond{i}, inputCondPath));
% % % %                 % Set the comment of the output file
% % % %                 sProcess.options.Comment.Value = sInputs(iInputCond(1)).Condition;
                % Process the average of condition #i
                tmpOutFiles = Concatenate(sProcess, sInputs(iInputCond));
                % Save the results
                OutputFiles = cat(2, OutputFiles, tmpOutFiles);
            end
        
        % === BY TRIAL GROUPS ===
        % Process each subject+condition+trial group independently
        case {5,6}
            % Get the condition path (SubjectName/Condition/CommentBase) or (Condition/CommentBase) for each input file
            CondPath = cell(1, length(sInputs));
            trialComment = cell(1, length(sInputs));
            for iInput = 1:length(sInputs)
                % Default comment
                trialComment{iInput} = sInputs(iInput).Comment;
                % If results/timefreq and attached to a data file
                if any(strcmpi(sInputs(iInput).FileType, {'results','timefreq'})) && ~isempty(sInputs(iInput).DataFile)
                    switch (file_gettype(sInputs(iInput).DataFile))
                        case 'data'
                            [sStudyAssoc, iStudyAssoc, iFileAssoc] = bst_get('DataFile', sInputs(iInput).DataFile);
                            if ~isempty(sStudyAssoc)
                                trialComment{iInput} = sStudyAssoc.Data(iFileAssoc).Comment;
                            else
                                bst_report('Warning', sProcess, sInputs(iInput), ['File skipped, the parent node has been deleted:' 10 sInputs(iInput).DataFile]);
                            end
                            
                        case {'results', 'link'}
                            [sStudyAssoc, iStudyAssoc, iFileAssoc] = bst_get('ResultsFile', sInputs(iInput).DataFile);
                            if ~isempty(sStudyAssoc)
                                [sStudyAssoc2, iStudyAssoc2, iFileAssoc2] = bst_get('DataFile', sStudyAssoc.Result(iFileAssoc).DataFile);
                                if ~isempty(sStudyAssoc2)
                                    trialComment{iInput} = sStudyAssoc2.Data(iFileAssoc2).Comment;
                                else
                                    bst_report('Warning', sProcess, sInputs(iInput), ['File skipped, the parent node has been deleted:' 10 sStudyAssoc.Result(iFileAssoc).DataFiles]);
                                end
                            else
                                bst_report('Warning', sProcess, sInputs(iInput), ['File skipped, the parent node has been deleted:' 10 sInputs(iInput).DataFile]);
                            end
                    end
                end
                % Subject average
                if (sProcess.options.avgtype.Value == 5)
                    CondPath{iInput} = [sInputs(iInput).SubjectName, '/', str_remove_parenth(trialComment{iInput})];
                % Grand average
                else
                    CondPath{iInput} = str_remove_parenth(trialComment{iInput});
                end
            end            
            uniquePath = setdiff(unique(CondPath), {''});
            % Process each condition path
            for i = 1:length(uniquePath)
                % Set progress bar at the same level for each loop
                if ~isempty(curProgress)
                    bst_progress('set', curProgress);
                end
                % Get all the files for condition #i
                iInputCond = find(strcmpi(uniquePath{i}, CondPath));
                % Do not process if there is only one input
                if (length(iInputCond) == 1)
                    bst_report('Warning', sProcess, sInputs(iInputCond(1)).FileName, 'File is alone in its trial/comment group. Not processed.');
                    continue;
                end
% % % %                 % Set the comment of the output file
% % % %                 sProcess.options.Comment.Value = str_remove_parenth(trialComment{iInputCond(1)});
                % Process the average of condition #i
                tmpOutFiles = Concatenate(sProcess, sInputs(iInputCond));
                % Save the results
                OutputFiles = cat(2, OutputFiles, tmpOutFiles);
            end
    end

end

function OutputFiles = Concatenate(sProcess, sInputs)

if strcmpi(sInputs(1).FileType,'data')
    OutputFiles = process_concat('Run',sProcess, sInputs);
    return;
end

isConnectn = ~isempty(strfind(sInputs(1).FileName,'connectn'));

% Load the first file, as the reference
NewMat = in_bst_timefreq(sInputs(1).FileName);
% Initialize time boundaries
StartTime = NewMat.Time(1);
EndTime   = NewMat.Time(end);
% Check all the input files
Ni = length(sInputs);
bst_progress('start', 'Processing', 'Concatenating...', 0, Ni);
for iInput = 2:Ni
    bst_progress('set', iInput);
    % Load the next file
    TimefreqMat = in_bst_timefreq(sInputs(iInput).FileName);
    % Check consistency with the number of sensors
    if (size(NewMat.TF,1) ~= size(TimefreqMat.TF,1))
        bst_report('Error', sProcess, sInputs(iInput), ['This file has a different number of channels than the previous ones: "' sInputs(iInput).ChannelFile '".']);
        return;
    end
    % Concatenate the F matrices
    NewMat.TF = [NewMat.TF, TimefreqMat.TF];
    
    % Update time boundaries
    StartTime = min(StartTime,TimefreqMat.Time(1));
    EndTime   = max(EndTime,TimefreqMat.Time(end));
end

% Set final time vector
NewMat.Time = linspace(StartTime,EndTime,size(NewMat.TF,2));

% % % Set comment
% % if (NewMat.Time(end) > 2)
% %     timeComment = sprintf('(%1.2fs,%1.2fs)', NewMat.Time(1), NewMat.Time(end));
% % else
% %     timeComment = sprintf('(%dms,%dms)', round(1000 * NewMat.Time(1)), round(1000 * NewMat.Time(end)));
% % end
% % NewMat.Comment = [str_remove_parenth(NewMat.Comment), ' concat' timeComment];

[sStudy, iStudy, Comment] = bst_process('GetOutputStudy', sProcess, sInputs);

% NewMat.Comment = [NewMat.Comment, ' | Concat ',Comment];
NewMat.Comment = [Comment, ' | concat'];

% 
NewMat.DataFile = [];
% Get output filename
fBase = 'timefreq_concat';
if isConnectn, fBase = 'timefreq_connectn_concat'; end
OutputFiles{1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), fBase);
% Save file
bst_save(OutputFiles{1}, NewMat, 'v6');
% Register in database
db_add_data(iStudy, OutputFiles{1}, NewMat);

end

