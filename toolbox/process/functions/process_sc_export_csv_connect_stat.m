function varargout = process_sc_export_csv_connect_stat( varargin )
% process_sc_export_csv_connect_stat: Export connectivity statistics in CSV
% file. We use 'import' inputtype since for now (oct-2013), you cannot
% slide statistic data file in process tab. This is why you have to enter
% the 'comment' of the file to compute (the name in the interface tree).
%  
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'EXPORT > CSV > Connectivity Statistics';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1765);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    sProcess.isSeparator = 1;
    % 
    sProcess.options.comment.Comment = 'Comment of inter-subject study to export:';
    sProcess.options.comment.Type    = 'text';
    sProcess.options.comment.Value   = [];
    % File selection options
    SelectOptions = {...
        '', ...                            % Filename
        '', ...                            % FileFormat
        'save', ...                        % Dialog type: {open,save}
        'Import BrainSuite folder...', ... % Window title
        'ExportData', ...                  % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                      % Selection mode: {single,multiple}
        'files', ...                        % Selection mode: {files,dirs,files_and_dirs}
        {'.csv', 'CSV (*.csv)', 'CSV'}, ... % Available file formats
        []};                          % DefaultFormats: {ChannelIn,DataIn,DipolesIn,EventsIn,MriIn,NoiseCovIn,ResultsIn,SspIn,SurfaceIn,TimefreqIn}
    % 
    sProcess.options.csvfile.Comment = 'CSV Exported file name:';
    sProcess.options.csvfile.Type    = 'filename';
    sProcess.options.csvfile.Value   = SelectOptions;
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

% Make sure it is connectivity statistic data
if isempty(strfind(sStat.FileName,'ptimefreq_'))||isempty(strfind(sStat.FileName,'connect'))
    strMsg = 'Inter-subject study must be connectivity statistics';
    bst_report('Error', sProcess, sInputs, strMsg);
    OutputFiles = {};
    return;
end

% Full name of CSV file to create
csvfile = sProcess.options.csvfile.Value{1};

% Get data struct
sMat = in_bst_timefreq(sStat.FileName);

% Generate cell table
pmax = 1;
[table, format] = ConnectivityStatTable(sMat,pmax);

% Write CSV file
WriteCSV(csvfile, table)

% === EXIT ===
OutputFiles = {sInputs.FileName}; % Return input file names

end

% ------------------------------------------------------------------------
function varargout = ConnectivityStatTable(s,pmax)
% Create a cell table of connectivity statistics values.
% 
%       ConnectivityStatTable(s,pmax)
%       T = ConnectivityStatTable(s,pmax)
%       [T,F] = ConnectivityStatTable(s,pmax)
% 
% T is a cell array table. First row is the header (columns names). Calling
% the function without output arguments prints the table in the command
% window. F is a cell array, whose elements contain the format
% (see FPRINTF) of each column of T.

if ~isfield(s,'RowNames') || ~isfield(s,'pmap')
    error('Input must be connectivity statistics');
end
if numel(s.RowNames)~=numel(s.RefRowNames)
    error('Unsupported 1xN connectivity yet');
end

if nargin<2, pmax = 1; end

[Nx,Nt,Nf] = size(s.pmap);

labelPairs = sc_bst_connect_vector_label(s.RowNames);
labelPairs = repmat(labelPairs(:),Nf,1);

Np = numel(labelPairs);

freqBand = cell(1,1,Nf);
freqBand(:) = s.Freqs(:,2);
freqBand = repmat(freqBand,[Nx,1,1]);
freqBand = freqBand(:);

pmap = s.pmap(:);
tmap = s.tmap(:);
df = s.df(:); if numel(df)==1, df = repmat(df,[Np,1]); end
mask = s.Options.statcond.mask(:);
mask_mp = s.Options.statcond.mask_mp(:);
alpha = s.Options.statcond.alpha;
mp = s.Options.statcond.mpmethod;

format   = {'%s',      '%.10f',       '%.10f',       '%d',                   '%d',                           '%.3f',      '%s'};
tablehdr = {'channels','p-values',    'score',       sprintf('p<%.3f',alpha),sprintf('p<%.3f (%s)',alpha,mp),'DoF',       'Band (Hz)'};
table    = [labelPairs,num2cell(pmap),num2cell(tmap),num2cell(mask),         num2cell(mask_mp),              num2cell(df),freqBand];

table = table(pmap<=pmax,:);

[sortval,sortidx] = sort(cell2mat(table(:,2)),'ascend');
table = table(sortidx,:);

% Final table
T = [tablehdr;table];

if nargout==0, disp(s.Comment); disp(T); end
if nargout>=1, varargout{1} = T; end
if nargout>=2, varargout{2} = format; end

end

% ------------------------------------------------------------------------
function WriteCSV(filename,M)

delimiter = ';';

% Make sure it is a cell matrix. The first row must be the header (strings)
if iscell(M)
    hdr = M(1,:);
% elseif isstruct(M)
%     hdrNames = fieldnames(M);
else
%     error('Input must be cell matrix or struct');
    error('Input must be cell matrix');
end

[Nr,Nc] = size(M);

% Determine the format (string, doubles, logical, ...) of each column
format = cell(1,Nc);
for iCol = 1:Nc
    Mi = M(2:end,iCol);
    if all(cellfun(@ischar, Mi))
        format{iCol} = '%s';
    elseif all(cellfun(@islogical, Mi)) || all(cellfun(@isinteger, Mi))
        format{iCol} = '%d';
    elseif all(cellfun(@isnumeric, Mi))
        format{iCol} = '%f';
    else
        error('Each column of cell matrix must be of unique type');
    end
end

% Open the file to write in it
fid = fopen(filename,'wb');

% Write the HEADER
iColFormat = ['%s',delimiter];
for iCol=1:Nc
    if iCol==Nc, iColFormat = ['%s']; end
    fprintf(fid, iColFormat, hdr{iCol});
end
fprintf(fid, '\r\n');

% Write the CSV data
for iRow = 2:Nr % Start at 2 since 1 is header
    for iCol = 1:Nc
        if iCol==Nc, iColFormat = [format{iCol}];
        else         iColFormat = [format{iCol},delimiter];
        end        
        fprintf(fid,iColFormat, M{iRow,iCol});
    end
    fprintf(fid, '\r\n');
end

fclose(fid);

end











