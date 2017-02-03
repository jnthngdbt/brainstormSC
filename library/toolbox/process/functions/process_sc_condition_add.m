function varargout = process_sc_condition_add( varargin )
% process_sc_condition_add: Add data to a condition
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'CONDITION > Add files';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(105);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    % 
    sProcess.options.condition.Comment = 'Condition name: ';
    sProcess.options.condition.Type    = 'text';
    sProcess.options.condition.Value   = 'NEWCOND';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>


% === INPUT PARAMETERS ===
Condition = sProcess.options.condition.Value;

% % Create condition for all subjects
% iStudy = db_add_condition('*', Condition, 1);
% if isempty(iStudy)
%     bst_report('Error', sProcess, sInputs, ['Cannot create condition : "' Condition '"']);
%     return;
% end

% === COMPUTE ===

Ni = numel(sInputs);

% Initialize progress bar
if bst_progress('isVisible'), bst_progress('set', bst_progress('get')); end
bst_progress('start', 'Processing', 'Adding files to condition...', 0, Ni);

% Add files to condition
OutputFiles = cell(Ni,1);
for iFile = 1:Ni
    bst_progress('set', iFile);

    FileMat = in_bst_data(sInputs(iFile).FileName);
    
    
    % ===== IMPORT FILES =====
    SubjectName = sInputs(iFile).SubjectName;
    % Get subject
    [sSubject, iSubject] = bst_get('Subject', SubjectName);
    % Define output study
    % Get condition asked by user
    [sStudy, iStudy] = bst_get('StudyWithCondition', fullfile(SubjectName, Condition));
    if isempty(sStudy)
        % Create condition 
        iStudy = db_add_condition(SubjectName, Condition, 1);
        if isempty(iStudy)
            bst_report('Error', sProcess, sInputs, ['Cannot create condition : "' Condition '"']);
            return;
        end
    end
    
    % Output filename
    OutputFile = strrep(file_fullpath(sInputs(iFile).FileName), sInputs(iFile).Condition, Condition);
    OutputFile = file_unique(OutputFile);
    % Save file
    bst_save(OutputFile, FileMat, 'v6');
    % Add file to database structure
    db_add_data(iStudy, OutputFile, FileMat);
    OutputFiles{iFile} = OutputFile;

    
end


end

