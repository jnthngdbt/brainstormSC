function varargout = process_sc_stat_convert2tf( varargin )
% process_sc_stat_convert2tf: Converts "statistics" typed files to 
% "timefrequency" typed files. Useful until Brainstorm support connectivity
% statistics. The field TF contains statistics values (t-values or 
% p-values). Newly created files are add in database in same study has
% input files. Unfriendlyly, Brainstorm cannot take statistics data in its
% process tab, so you have to manually enter the name of the file to
% process.
% 
% External call:
% 
%   OutputFiles = process_sc_stat_convert2tf('Stat2TF',sInputs)
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'STATISTICS > CONNECTIVITY > Convert "stat" to "time-freq" (for visualization)';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1555);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    sProcess.isSeparator = 1;
    % 
    sProcess.options.comment.Comment = 'Enter the name of inter-subject study to convert:';
    sProcess.options.comment.Type    = 'text';
    sProcess.options.comment.Value   = [];
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

% Get inter-subject studies
[sStudy, iStudy] = bst_get('AnalysisInterStudy');

% Find the study to export, selected by the user using its 'comment'
comment = sProcess.options.comment.Value;
iMatch = find(strcmpi(comment,{sStudy.Stat.Comment}));
if isempty(iMatch)
    strMsg = sprintf('No inter-subject study has comment [%s]',comment);
    bst_report('Error', sProcess, sInputs, strMsg);
    OutputFiles = {};
    return;
elseif numel(iMatch)>1
    strMsg = sprintf('Inter-subject study comment [%s] is not unique. Taking first study.',comment);
    bst_report('Error', sProcess, sInputs, strMsg);
    iMatch = iMatch(1);
end
sStat = sStudy.Stat(iMatch);

sStat.iStudy = iStudy;

% Make sure it is connectivity statistic data
if isempty(strfind(sStat.FileName,'ptimefreq_'))||isempty(strfind(sStat.FileName,'connect'))
    strMsg = 'Inter-subject study must be connectivity statistics';
    bst_report('Error', sProcess, sInputs, strMsg);
    OutputFiles = {};
    return;
end

OutputFiles = Stat2TF( sStat );

end

function OutputFiles = Stat2TF( sInputs )
% Stat2TF: Converts "statistics" typed files to 
% "timefrequency" typed files. Useful until Brainstorm support connectivity
% statistics. The field TF contains statistics values (t-values or 
% p-values). Newly created files are add in database in same study has
% input files.
%
% Authors: Jonathan Godbout, 2013

if isfield(sInputs,'FileType')
    if ~strcmpi(unique({sInputs.FileType}),'ptimefreq')
        error('All files must be timefrequency statistics files');
    end
end

OutputFiles = cell(numel(sInputs),1);
for iInput = 1:numel(sInputs)
    
    % Load statistics data
    sStat = in_bst_timefreq(sInputs(iInput).FileName);
    % Create time-frequency data from stat
    sTF = db_template('timefreqmat');
    for iField=intersect(fieldnames(sStat),fieldnames(sTF))'
        sTF.(iField{1}) = sStat.(iField{1});
    end
    sTF.TF          = sStat.tmap;
    sTF.DataType    = 'data';
    sTF.Options.statconvert = sStat;
%     sTF.Measure = 'other';
    sTF.Method = 'cohere';  % Otherwise display bug
    
    isConnectn = ~isempty(sTF.RefRowNames);
    if isConnectn,  FileType = 'timefreq_connectn';
    else            FileType = 'timefreq';
    end
    
    % ===== SAVE FILE =====
    % Study
    sStudy = bst_get('Study', sInputs(iInput).iStudy);
    % Output filename
    OutputFiles{iInput} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), FileType);
    % Save on disk
    bst_save(OutputFiles{iInput}, sTF, 'v6');
    % Register in database
    db_add_data(sInputs(iInput).iStudy, OutputFiles{iInput}, sTF);
end
    
% Update database interface
db_reload_studies([sInputs.iStudy]);

end
