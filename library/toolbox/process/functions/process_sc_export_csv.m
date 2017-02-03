function varargout = process_sc_export_csv( varargin )
% process_sc_export_csv: Export data to CSV. Multiple variables (condition,
% subject, channels (pairs for connectivity), values...). It only creates 1
% CSV file; multiple inputs will be merged (should all be of same data
% type).
%  
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'EXPORT > CSV';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1760);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data','timefreq'};
    sProcess.OutputTypes = {'data','timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    % Specific sensor pairs
    sProcess.options.pairfilter.Comment = 'Only export specific sensors pairs (ex: F3-P3, F3-F4). Connectivity data only.';
    sProcess.options.pairfilter.Type    = 'text';
    sProcess.options.pairfilter.Value   = [];
    % Specific frequency bands
    sProcess.options.freqfilter.Comment = 'Only export specific frequency bands (All=0)';
    sProcess.options.freqfilter.Type    = 'value';
    sProcess.options.freqfilter.Value   = {[1 2],'',0};
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

% Get current progressbar position
if bst_progress('isVisible')
    curProgress = bst_progress('get');
else
    curProgress = [];
end

% Make sure it is all the same data type
FileType = unique({sInputs.FileType});
if numel(FileType)>1
    strMsg = sprintf('All data inputs should be of same type');
    bst_report('Error', sProcess, sInputs, strMsg);
    OutputFiles = {};
    return;
else
    FileType = FileType{1};
    % Special 'timefreq' case: connectivity
    switch lower(FileType)
        case 'timefreq'
            sMat = in_bst_timefreq(sInputs(1).FileName);
            if ~isempty(sMat.RefRowNames)
                FileType = 'connect';
            end
    end
end

% Full name of CSV file to create
csvfile = sProcess.options.csvfile.Value{1};

% Built a table matrix from all files
bst_progress('text', 'Exporting files to CSV...');
Ni = numel(sInputs);
for iInput=1:Ni
    bst_progress('set',  round(curProgress + 100*iInput/Ni));
    
    SubjectName = sInputs(iInput).SubjectName;
    Condition = sInputs(iInput).Condition;
    
    % Create cell matrix table of the data
    switch lower(FileType)
        case 'data'
            error('Filteype [%s] support TODO', FileType)
        case 'timefreq'
            error('Filteype [%s] support TODO', FileType)
        case 'connect'
            % Get data struct
            sMat = in_bst_timefreq(sInputs(iInput).FileName);
            if ~isreal(sMat.TF)
                bst_report('Error',sProcess, sInputs(iInput),'Data is complex! Must apply a measure before exporting.');
                continue;
            end
            Ti = ConnectivityTable(sMat, sProcess.options.freqfilter.Value{1}, sProcess.options.pairfilter.Value);
        otherwise
            error('Filteype [%s] not supported', FileType)
    end
    
    % Deal with header
    if iInput==1 % Create header
        T = Ti(1,:);
        T = [{'subject'}, {'condition'}, T];
    end 
    Ti(1,:) = []; % Remove header 
    
    % Table's number of rows and columns
    [Nr, Nc] = size(Ti);
    
    % Add subject and condition columns
    Ti = [repmat({SubjectName},Nr,1), repmat({Condition},Nr,1), Ti];
    
    % Add to master table 
    T = [T; Ti];
    
end

% Write CSV file
WriteCSV(csvfile, T);

% === EXIT ===
OutputFiles = {sInputs.FileName}; % Return input file names

end

% ------------------------------------------------------------------------
function varargout = ConnectivityTable(sMat, iFreqFilter, PairFilter)
% Create a cell table of connectivity values.
% 
%       ConnectivityStatTable(sMat)
%       T = ConnectivityStatTable(sMat)
%       [T,F] = ConnectivityStatTable(sMat)
% 
% T is a cell array table. First row is the header (columns names). Calling
% the function without output arguments prints the table in the command
% window. F is a cell array, whose elements contain the format
% (see FPRINTF) of each column of T.

if nargin<2, iFreqFilter = 0; end
if nargin<3, PairFilter = []; end

% See if data is supported
if ~isfield(sMat,'RowNames')
    error('Input must be connectivity');
else
    if isempty(sMat.RefRowNames)
        error('Input must be connectivity');
    end
end
if numel(sMat.RowNames)~=numel(sMat.RefRowNames)
    error('Unsupported 1xN connectivity yet');
end
[Nx,Nt,Nf] = size(sMat.TF);
Nc = numel(sMat.RowNames);
if Nt>1
    error('Unsupported dynamic connectivity yet (with time dimension)');
end
if Nx~=(Nc*(Nc+1)/2)
    error('Unsupported directed connectivity yet');
end

% Frequency band filter
if iFreqFilter==0, iFreqFilter = 1:size(sMat.TF,3); end

% Channel pairs filter
AllPairs = sc_bst_connect_vector_label(sMat.RowNames);
if isempty(PairFilter)
    PairFilter = AllPairs;
else
    PairFilter = textscan(PairFilter,'%s','delimiter',',');
    PairFilter = PairFilter{1};
end
PairMatrix = sc_bst_connect_matrix_label(sMat.RowNames);
iPairFilter = false(size(PairMatrix));
for i=1:numel(PairFilter)
    iPairFilter = iPairFilter | strcmpi(PairMatrix,PairFilter{i});
end
iPairFilter = iPairFilter | iPairFilter'; % Undirected connectivity
iPairFilter = find(sc_bst_connect_format_mat2vec(iPairFilter)>0);

% Filter TF
TF = sMat.TF(iPairFilter,:,iFreqFilter);

[Nx,Nt,Nf] = size(TF);

% Number of table rows (without header)
Nr = Nx*Nt*Nf; 

% Strings of channel pairs
% labelPairs = sc_bst_connect_vector_label(sMat.RowNames);
labelPairs = AllPairs(iPairFilter);
labelPairs = repmat(labelPairs(:),Nf,1);

% Frequencies
if iscell(sMat.Freqs)
    Freqs = cell(1,1,Nf);
    Freqs(:) = sMat.Freqs(iFreqFilter,2);
    Freqs = repmat(Freqs,[Nx,1,1]);
    Freqs = Freqs(:);
else
    Freqs = cell(1,1,Nf);
    Freqs(:) = cellstr(num2str(sMat.Freqs(iFreqFilter)));
    Freqs = repmat(Freqs,[Nx,1,1]);
    Freqs = Freqs(:);
end

% Method
Method = repmat({sMat.Method},Nr,1);
% Comment
Comment = repmat({sMat.Comment},Nr,1);

% % % % % format   = {'%s',       '%s',            '%s'     '%s'     '%.10f'};
% % % % % tablehdr = {'channels', 'frequency (Hz)' 'method' 'comment' 'values'};
% % % % % table    = [labelPairs, Freqs,           Method,  Comment,  num2cell(sMat.TF(:))];

format   = {'%s',       '%s',            '%s'     '%.10f'};
tablehdr = {'channels', 'frequency (Hz)' 'comment' 'values'};
table    = [labelPairs, Freqs,           Comment,  num2cell(TF(:))];

% Final table
T = [tablehdr;table];

if nargout==0, disp(sMat.Comment); disp(T); end
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











