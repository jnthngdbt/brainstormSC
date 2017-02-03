function varargout = process_sc_import_harmonie( varargin )
% process_sc_import_harmonie: Convert Stellate Harmonie (*.sts, *.sig)
% dataset into a Brainproducts dataset (*.vhdr, *.vmrk, *.eeg) that 
% Brainstorm can handle. Then imports the latter with built-in BST import
% functions.
% 
% External function calls to convert Stellate->Brainproducts
% 
%       process_sc_import_harmonie('Harmonie2Brainproducts')
%       process_sc_import_harmonie('Harmonie2Brainproducts',FILES);
%       process_sc_import_harmonie('Harmonie2Brainproducts',..., 'addchannelinfo');
%
% Author: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'CONVERT-IMPORT > Stellate Harmonie (only under Windows and Matlab 32bits)';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(16);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    sProcess.isSeparator = 0;
    % Option: File to import
    sProcess.options.datafile.Comment = 'You will be asked to select Stellate files';
    sProcess.options.datafile.Type    = 'label';
    sProcess.options.datafile.Value   = [];
    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
    % Channel
    sProcess.options.addchaninfo.Comment = 'Add channel info in event names';
    sProcess.options.addchaninfo.Type    = 'checkbox';
    sProcess.options.addchaninfo.Value   = 0;
    % Align sensors
    sProcess.options.channelalign.Comment = 'Align sensors using headpoints';
    sProcess.options.channelalign.Type    = 'checkbox';
    sProcess.options.channelalign.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

datafiles = sProcess.options.datafile.Value;

if isempty(datafiles)
    % Get file names
    [datafiles, pname] = uigetfile( ...
    {'*.sig;*.SIG' , 'Stellate Harmonie Files (*.sig)'}, ...
    'Pick Stellate Harmonie files (YOU MUST BE UNDER WINDOWS AND MATLAB 32BITS)', ...
    'MultiSelect', 'on');
    if ~iscell(datafiles), datafiles = {datafiles}; end
    for i=1:numel(datafiles)
        datafiles{i} = [pname,datafiles{i}];
    end
end
if ~iscell(datafiles), datafiles = {datafiles}; end
Nf = numel(datafiles);

% Loop over files
OutputFiles = cell(Nf,1);
for iInput=1:Nf

    % File name infos
    [p,subjectname,ext] = fileparts(datafiles{iInput});
    
    OutputFiles{iInput} = [datafiles{iInput},'.eeg'];
    
    try
        % Skip conversion if already done (converted file exists)
        if exist(OutputFiles{iInput},'file')>0
            strMsg = sprintf('Stellate file [%s] is already converted to Brainproduct.',[subjectname,ext]);
            bst_report('Warning', sProcess, sInputs, strMsg)
        else
            if sProcess.options.addchaninfo.Value>0 % Add channel in event name
                Harmonie2Brainproducts(datafiles{iInput},'addchannelinfo');      
            else
                Harmonie2Brainproducts(datafiles{iInput});      
            end
        end
        % Import in BST database
        % Process: Create link to raw file
        sFiles = bst_process(...
            'CallProcess', 'process_import_data_raw', ...
            [], [], ...
            'subjectname', subjectname, ...
            'datafile', {OutputFiles{iInput}, 'EEG-BRAINAMP'}, ...
            'channelalign', sProcess.options.channelalign.Value);
        % Overwrite with content of sFiles, otherwise file is not recognized
        OutputFiles{iInput} = sFiles.FileName;

    catch ME
        OutputFiles{iInput} = [];
        strMsg = sprintf('There was a problem converting file [%s]: \n %s',[subjectname,ext],getReport(ME));
        bst_report('Warning', sProcess, sInputs, strMsg)
    end
end
% % Remove empty output files
% OutputFiles(cellfun(@isempty,OutputFiles)) = [];

end

% ========================================================================
% ========================================================================
function Harmonie2Brainproducts(varargin)
% Convert Stellate files (*.sts, *.sig) to Brainproduct format (*.eeg,
% *.vhdr, *.vmrk).
% 
%   Harmonie2Brainproducts;
%   Harmonie2Brainproducts(FILES);
%   Harmonie2Brainproducts(..., 'addchannelinfo');
% 
% where FILES is a CELL array of strings (Stellate file names) or a string
% (single Stellate file). With no input, will ask user to select the files
% to convert. The 'addchannelinfo' flag will add the channel name of the
% markers in their "description" field (only if the marker is on a specific
% channel; "description" wont be changed if the marker is on all channels).

% Option flag: add channel info in name
iAddChannel = find(strcmpi(varargin,'addchannelinfo'));
if ~isempty(iAddChannel)
    varargin(iAddChannel) = [];
    iAddChannel = true;
else
    iAddChannel = false;
end

[fname,pname] = get_files(varargin{:});
if isempty(fname), return; end

Nf = numel(fname);

for ii=1:Nf

    [hdr,mrk] = OpenHarmonie([pname,fname{ii}]);
    hdr = write_vhdr(pname, fname{ii}, hdr);
    mrk = write_vmrk(pname, fname{ii}, hdr, mrk, iAddChannel);

    write_eeg(pname, fname{ii}, hdr);
    
end

end

% =========================================================================
function [fname,pname] = get_files(varargin)

if ~isempty(varargin)
    if isempty(varargin{1}), varargin = []; end
end

fname = []; % Filenames (1xN CELL)
pname = []; % Pathname (CHAR)
if isempty(varargin)
    [fname, pname] = uigetfile( ...
                     {'*.sig;*.SIG' , 'Stellate Harmonie SIG files (*.sig)'}, ...
                      'Pick the SIG files to convert', ...
                      'MultiSelect', 'on');

    if ~iscell(fname), fname = {fname}; end
elseif iscell(varargin{1})
    %%%%%%
    fname = varargin{1};
    for i=1:numel(fname)
        [p,n,e] = fileparts(fname{i}); if ~isempty(p), p = [p,filesep]; end
        pname = p;
        fname{i} = [n,e];
    end
elseif ischar(varargin{1})
    %%%%%%%
    [p,n,e] = fileparts(varargin{1}); if ~isempty(p), p = [p,filesep]; end
    pname = p;
    fname = {[n,e]};
else
    error('I dont know what you want to do');
end

end
% =========================================================================
function [hdr mrk] = OpenHarmonie(varargin)

[fname,pname] = get_files(varargin);
if isempty(fname), return; end
fname = fname{1};

hdr = ReadHarmonieHeader([pname,fname]);
mrk = ReadHarmonieMarker([pname,fname]);

hdr.date = get_date_stellate(hdr,mrk);

for i=1:numel(mrk)
    mrk(i).feature = add_features(hdr,mrk(i).feature);
    % Detects different file labels and channel labels (differential montage)
    if all(mrk(i).feature.channel_i)==0
        for k=1:numel(mrk(i).feature.channel_i)
            idx = strfind(mrk(i).feature.channel{k},'-');
            if ~isempty(idx)
                mrk(i).feature.channel{k} = mrk(i).feature.channel{k}(1:idx-1);
            end
        end
        mrk(i).feature = add_features(hdr,mrk(i).feature);
    end
end

end

function mrk = add_features(hdr,mrk)

% Add CHANNEL_I feature from channel labels
label = hdr.label;
N = length(mrk.channel);
mrk.channel_i = zeros(1,N);
for i=1:N
    idx = find(strcmpi(mrk.channel{i},label));
    if ~isempty(idx)
        mrk.channel_i(i) = idx(1);
    end
end

end

% -------------------------------------------------------------------------
function date = get_date_stellate(hdr,mrk)

date_num_off = datenum(hdr.orig.RecordingStartTime,'yyyy-mm-dd HH:MM:SS');
offset = get_offset(mrk);
secs_in_day = 24*60*60;
date_offset = (offset/hdr.Fs)/secs_in_day;
date_num = date_num_off - date_offset;

date.vec = datevec(date_num);
date.str = datestr(date_num,'yyyy-mm-dd HH:MM:SS.FFF');
date.num = date_num;

end

function offset = get_offset(mrk)

offset = 0;
for ii=1:numel(mrk)
    switch mrk(ii).name
        case 'Section d''échantillon.'
            if mrk(ii).feature.position_i(1)>0
                offset = mrk(ii).feature.position_i(1)-1;
            end
        case 'Sample Section'
            if mrk(ii).feature.position_i(1)>0
                offset = mrk(ii).feature.position_i(1)-1;
            end
        otherwise
    end
end

end

% =========================================================================
function hdr = ReadHarmonieHeader(varargin)

% if ~strcmpi('PCWIN',computer)
%     error('Handling Stellate files (.sig, .sts) is only available on Windows 32 bits.')
% end

[fname,pname] = get_files(varargin);
fname = fname{1};
Filename = [pname,fname];

if ~exist('mFileOpen','file')
    error('Stellate-Matlab interface not in path!');
end

% Try to open SIG file using DLL interface
try
    mFileOpen(Filename);
catch
    error('Handling Stellate files (.sig, .sts) is only available on Windows and Matlab 32 bits.')
%     error('Could not use Stellate-Matlab interface DLL file. This can be an operating system problem (Win32, Win64, OSX, ...) or a Matlab version problem (may only work on 2009-2010 and older). ');
end

% check the file length
[orig.NumRecs, orig.NumSamps, orig.NumSecs] = mGetFileLength(Filename);

% check the number of recording channels
[orig.NumOfChan] = mGetNumChan(Filename);

% check the base sampling frequency
[orig.TrueSampFreq] = mGetTrueSampFreq(Filename);

% the time of the first recorded sample in the given SIG file.
[orig.RecordingStartTime] = mGetRecStartTime(Filename);

% read the montage list
[orig.MtgList] = mGetMtgList(Filename);

% read the electrode names of the recording montage
[orig.RecChanNames] = mGetRecChanName(Filename);

% check the sampling rate of every recording channel 
[orig.SampRate] = mGetChanSampRate(Filename);

% check the event group list
[orig.GrpList] = mGetEvtGrpList(Filename);

mFileClose;

%              Fs: 250
%          nChans: 67
%           label: {67x1 cell}
%        nSamples: 1245489
%     nSamplesPre: 0
%         nTrials: 1
%            orig: [1x1 struct]

hdr.Fs = round(orig.TrueSampFreq);
hdr.nChans = orig.NumOfChan;
hdr.label = orig.RecChanNames;
hdr.nSamples = orig.NumSamps - 1; % Crash when accessing last sample...
hdr.nSamplesPre = 0;
hdr.nTrials = orig.NumRecs;
hdr.orig = orig;
end

% =========================================================================
function evt = ReadHarmonieMarker(varargin)
% A lot of limitations. Only on Windows.

% if ~strcmpi('PCWIN',computer)
%     error('Handling Stellate files (.sig, .sts) is only available on Windows 32 bits.')
% end

[fname,pname] = get_files(varargin);
fname = fname{1};

Filename = [pname,fname];

if ~exist('mFileOpen','file')
    error('Stellate-Matlab interface not in path!');
end

% Try to open SIG file using DLL interface
try
    mFileOpen(Filename);
catch
    error('Handling Stellate files (.sig, .sts) is only available on Windows and Matlab 32 bits.')
%     error('Could not use Stellate-Matlab interface DLL file. This can be an operating system problem (Win32, Win64, OSX, ...) or a Matlab version problem (may only work on 2009-2010 and older). ');
end

% check the base sampling frequency
fs = mGetTrueSampFreq(Filename);
fs = round(fs);

% check the event group list
[GrpList] = mGetEvtGrpList(Filename);

Ng = numel(GrpList);

mrk = cell(1,Ng);
mrk(:) = {struct};

for ii=1:numel(GrpList)
    % retrieves all the status items and their properties 
    [PropHeader, StatusItems, channel] = mGetStatusItems(Filename, GrpList{ii});
    
    if isempty(StatusItems)
        mrk{ii}.type = GrpList(ii);
        mrk{ii}.position_i = 0;
        mrk{ii}.duration_i = 0;
        mrk{ii}.channel = {'ALL'};
        continue;
    elseif isempty(channel{1}), channel(:) = {'ALL'};
    end
    
    [Nv Nf] = size(StatusItems); % Number of values (1) and fields (2)
    
    fields = textscan(PropHeader,'%s',Nf);
    fields = fields{:};
    
    N = mGetNumStatusItemsOfEvt(Filename,GrpList{ii});
    
    mrk{ii}.type = cell(1,Nv);
    mrk{ii}.type(:) = GrpList(ii);
    mrk{ii}.channel = channel(:)';
    for ff=1:Nf
        switch lower(fields{ff})
            case 'beginsample'
                mrk{ii}.position_i = StatusItems(:,ff)';
            case 'eventinterval'
                mrk{ii}.duration_i = StatusItems(:,ff)';
            case 'type'
                continue;
            case 'channel'
                continue;
            otherwise
                s = lower(fields{ff});
                s(strfind(s,'é')) = 'e';
                s(strfind(s,'è')) = 'e';
                s(strfind(s,'à')) = 'a';
                s(strfind(s,'/')) = '';
                s(strfind(s,'\')) = '';
                s(strfind(s,'-')) = '';
                s(strfind(s,'.')) = '';
                s(strfind(s,' ')) = '';
                s(strfind(s,':')) = '';
                if ismember(s(1),'0123456789'), s = ['x',s]; end
                if isempty(s), continue; end
                mrk{ii}.(s) = StatusItems(:,ff)';
        end
    end
    mrk{ii}.position_i(mrk{ii}.position_i==0) = 1;
end

mFileClose;

mrk = normalize_marker_stellate(mrk);

% Create an array of structures
evt = [];
for ii=1:numel(GrpList)
    evti.name = mrk{ii}.type{1};
    evti.feature = mrk{ii};
    evt = [evt,evti];
end

end

function mrk = normalize_marker_stellate(mrk)
% MRK is a CELL array of struct obtained from READ_MARKER_STELLATE

for i=1:numel(mrk)
    mrk{i}.description = mrk{i}.type; % Add description feature
    type = unique(mrk{i}.type);
    switch lower(type{1})
        case {'section d''échantillon.', 'sample section'}
            mrk{i}.type(:) = {'Sample Section'};
            mrk{i}.description = mrk{i}.type;
        case {'stage','stade'}
            mrk{i}.type(:) = {'Stage'};
            sleepDescription = getSleepDescription;
            mrk{i}.description = sleepDescription(mrk{i}.stage+1);
        case {'m-eve','m-eves','m-éve','m-éves','micro-éveil','micro-éveils'}
            mrk{i}.type(:) = {'m-eve'};
            mrk{i}.description = mrk{i}.type;
    end
    
    if ~isempty(strfind(lower(type{1}),'harmact'))
                mrk{i}.type(:) = {'harmact'};
                mrk{i}.description = mrk{i}.type;
    end
end

end

function sleepDescription = getSleepDescription
sleepDescription = {'WAKE','NREM1','NREM2','NREM3','NREM4','REM','NA','NA','NA','NA'};
end

% =========================================================================
function hdr = write_vhdr(p, n, hdr)

% hdr.orig = read_header_stellate([p,n]);
hdr.orig = hdr;

hdr.Codepage            = 'UTF-8';
hdr.DataFile            = [n '.eeg'];
hdr.MarkerFile          = [n '.vmrk'];
hdr.DataFormat          = 'BINARY';
hdr.DataOrientation     = 'MULTIPLEXED';
hdr.DataType            = 'TIMEDOMAIN';
hdr.NumberOfChannels  	= hdr.orig.nChans;
hdr.DataPoints          = hdr.orig.nSamples - 1; %%%%%% Accessing last point crashes... SEE ALSO IN WRITE EEG
hdr.SamplingInterval    = 10^6/hdr.orig.Fs;
hdr.BinaryFormat        = 'IEEE_FLOAT_32';
hdr.UseBigEndianOrder   = 'NO';
hdr.resolution          = ones(hdr.orig.nChans,1);
hdr.reference           = cell(hdr.orig.nChans,1); %***
hdr.label               = hdr.orig.label;
hdr.units               = [char(181) 'V'];

% open the header file and write the ascii header information
fid = fopen([p n '.vhdr'], 'wb');
fprintf(fid, 'Brain Vision Data Exchange Header File Version 1.0\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '[Common Infos]\r\n');
fprintf(fid, 'Codepage=%s\r\n',          hdr.Codepage);
fprintf(fid, 'DataFile=%s\r\n',          hdr.DataFile);
fprintf(fid, 'MarkerFile=%s\r\n',        hdr.MarkerFile);
fprintf(fid, 'DataFormat=%s\r\n',        hdr.DataFormat);
fprintf(fid, 'DataOrientation=%s\r\n',   hdr.DataOrientation);
fprintf(fid, 'DataType=%s\r\n',          hdr.DataType);
fprintf(fid, 'NumberOfChannels=%d\r\n',  hdr.NumberOfChannels);
fprintf(fid, 'DataPoints=%d\r\n',        hdr.DataPoints);
fprintf(fid, 'SamplingInterval=%.15g\r\n',  hdr.SamplingInterval);
fprintf(fid, '\r\n');
fprintf(fid, '[Binary Infos]\r\n');
fprintf(fid, 'BinaryFormat=%s\r\n',      hdr.BinaryFormat);
fprintf(fid, 'UseBigEndianOrder=%s\r\n', hdr.UseBigEndianOrder);
fprintf(fid, '\r\n');
fprintf(fid, '[Channel Infos]\r\n');
fprintf(fid, '; Each entry: Ch<Channel number>=<Name>,<Reference channel name>,\r\n');
fprintf(fid, '; <Resolution in "Unit">,<Unit>, Future extensions...\r\n');
fprintf(fid, '; Fields are delimited by commas, some fields might be omited (empty).\r\n');
fprintf(fid, '; Commas in channel names are coded as "\1".\r\n');
for cc=1:hdr.NumberOfChannels
  fprintf(fid, 'Ch%d=%s,,%g,%s\r\n', cc, hdr.label{cc}, hdr.resolution(cc),hdr.units);
end

fclose(fid);

end

% =========================================================================
function mrk = write_vmrk(p,n,hdr,mrk,iAddChannel)

if nargin<5, iAddChannel = false; end

% mrk = read_marker_stellate([p,n],hdr);
mrk = convertMrk2brainproducts(mrk);

if iAddChannel % Add channel name in marker descriptions
    mrk = add_channel_info(mrk,hdr);
end

% Open the marker file and write the ascii header information
fid = fopen([p n '.vmrk'], 'wb');
fprintf(fid, 'Brain Vision Data Exchange Marker File, Version 1.0\r\n');
fprintf(fid, '\r\n');
fprintf(fid, '[Common Infos]\r\n');
fprintf(fid, 'Codepage=%s\r\n',     'UTF-8');
fprintf(fid, 'DataFile=%s\r\n',     [n '.eeg']);
fprintf(fid, '\r\n');

fprintf(fid, '[Marker Infos]\r\n');
fprintf(fid, '; Each entry: Mk<Marker number>=<Type>,<Description>,<Position in data points>,\r\n');
fprintf(fid, '; <Size in data points>, <Channel number (0 = marker is related to all channels)>\r\n');
fprintf(fid, '; Fields are delimited by commas, some fields might be omitted (empty).\r\n');
fprintf(fid, '; Commas in type or description text are coded as "\1".\r\n');

% The first marker is always a 'New Segment'. However, it is not present in
% the |.mat| file, so it has to be created. In its case, the field <Channel
% number> corresponds to the date.
% % % % date_str = get_date_stellate(hdr,mrk);
date_str = hdr.orig.date.str;
date_str(findstr(date_str,'.')) = [];
date_str(findstr(date_str,':')) = [];
date_str(findstr(date_str,'-')) = [];
date_str(findstr(date_str,'/')) = [];
date_str(findstr(date_str,'\')) = [];
date_str(findstr(date_str,' ')) = [];

% Write the 'New Segment' with the date string
fprintf(fid, 'Mk%d=%s,%s,%d,%d,%d,%s\r\n', 1, 'New Segment', ...
               '', 1, 1, 0, date_str);

try
    N = numel(mrk.type);
catch
    fclose(fid); return;
end

for mm=1:N
    if strcmpi(mrk.type{mm},'New Segment')
        n = fieldnames(mrk);
        for ii=1:numel(n), mrk.(n{ii})(mm) = []; end
        break;
    end
end

N = numel(mrk.type);
% For each marker, get its info and write in file
for mm=1:N
    fprintf(fid,'Mk%d=%s,%s,%d,%d,%d\r\n',mm+1,...
                                          mrk.type{mm}, ...
                                          mrk.description{mm}, ...
                                          mrk.start_i(mm), ...
                                          mrk.duration_i(mm), ...
                                          mrk.channel_i(mm));
end

fclose(fid);

end

function mrk = convertMrk2brainproducts(mrk)

mrkOut.type = [];
mrkOut.description = [];
mrkOut.start_i = [];
mrkOut.duration_i = [];
mrkOut.channel_i = [];
for i=1:numel(mrk)
    mrkOut.type         = [mrkOut.type, mrk(i).feature.type];
    mrkOut.description  = [mrkOut.description, mrk(i).feature.description];
    mrkOut.start_i      = [mrkOut.start_i, mrk(i).feature.position_i];
    mrkOut.duration_i   = [mrkOut.duration_i, mrk(i).feature.duration_i];
    mrkOut.channel_i    = [mrkOut.channel_i, mrk(i).feature.channel_i];
end

[sortv,sorti] = sort(mrkOut.start_i,'ascend');
mrkOut.type         = mrkOut.type(sorti);
mrkOut.description  = mrkOut.description(sorti);
mrkOut.start_i      = mrkOut.start_i(sorti);
mrkOut.duration_i   = mrkOut.duration_i(sorti);
mrkOut.channel_i    = mrkOut.channel_i(sorti);


mrk = mrkOut;


end

function mrk = add_channel_info(mrk,hdr)

% Markers not on all channels
iNotAll = find(mrk.channel_i~=0); 

% If all markers are on all channels, nothing to do: return
if isempty(iNotAll), return; end

% Create a duplicate set of markers with added channel info
mrknew = mrk;
for i=1:numel(iNotAll)
    mrknew.description{iNotAll(i)} = sprintf('%s_%s', ...
    mrknew.description{iNotAll(i)}, hdr.label{mrknew.channel_i(iNotAll(i))});
end

% Append newly created markers to original set
mrk.type         = [mrk.type,        mrknew.type(iNotAll)];
mrk.description  = [mrk.description, mrknew.description(iNotAll)];
mrk.start_i      = [mrk.start_i,     mrknew.start_i(iNotAll)];
mrk.duration_i   = [mrk.duration_i,  mrknew.duration_i(iNotAll)];
mrk.channel_i    = [mrk.channel_i,   mrknew.channel_i(iNotAll)];

% Sort all markers in chronological order
[sortv,sorti] = sort(mrk.start_i,'ascend');
mrk.type         = mrk.type(sorti);
mrk.description  = mrk.description(sorti);
mrk.start_i      = mrk.start_i(sorti);
mrk.duration_i   = mrk.duration_i(sorti);
mrk.channel_i    = mrk.channel_i(sorti);

end

% =========================================================================
function write_eeg(p,n,hdr)
% Accessing the last sample seems to crash...

length_i = min(2^16, hdr.DataPoints)-1;

filename = [p,n,'.eeg'];

fid = fopen(filename,'w');
N = floor(hdr.DataPoints / length_i);
start_i = 1;
[data, SampRate] = mGetDataSamp([p,n], start_i, length_i);
data = interpolateLowerSamprate(data,SampRate,hdr);
try bst_progress('start', 'Converting Stellate Harmonie -> Brainproducts', 'Writing data...', 0, N); end
for ii=1:N
    try bst_progress('set', ii); end
    fprintf('%.5f %%\n',100*ii/N);
   
    fwrite(fid,data','float32'); % Transpose data!
    
    start_i = start_i + length_i;
    if ii<N, 
        [data, SampRate] = mGetDataSamp([p,n], start_i, length_i);
        data = interpolateLowerSamprate(data,SampRate,hdr);
    end
end
% Write remaining samples
length_i = hdr.DataPoints - start_i + 1; % Accessing the last sample seems to crash... BUT ALREADY SUBTRACTED ONE POINT IN WRITING HEADER
if length_i>0
    [data, SampRate] = mGetDataSamp([p,n], start_i, length_i);
    data = interpolateLowerSamprate(data,SampRate,hdr);
    fwrite(fid,data','float32');% Transpose data!
end
fclose(fid);

end

function data = interpolateLowerSamprate(data,SampRate,hdr)

SampRate = round(SampRate);
% fs = hdr.orig.orig.TrueSampFreq;
fs = round(hdr.orig.orig.TrueSampFreq); %%%
toInterp = find(SampRate~=fs);
N = size(data,1);
Ni = floor(N*SampRate/fs);
for i=1:numel(toInterp)
    idx = toInterp(i);
    data(:,idx) = interpft(data(1:Ni(idx),idx),N);
end

end

