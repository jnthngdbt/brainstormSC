function varargout = process_sc_import_edf( varargin )
% process_sc_import_edf: Create link to raw EDF file. The difference with
% built-in BST EDF import (July-2013) is that it deals with channels having
% multiple sampling frequencies; if it is the case, it converts EDF files
% to other EDF files, interpolating channels with lower sampling
% frequencies to the highest sampling frequency, so that BST can import and
% use these new files.
% 
%   OutputFiles = process_sc_import_edf('Run', sProcess, sInputs)
%   hdrUSR = process_sc_import_edf('Convert2USR', hdrMSR)
%   hdr = process_sc_import_edf('ReadEDFHeader', filename)
%   F = process_sc_import_edf('ReadEDFData', hdr, TimeBounds, iChannels)
%
% Author: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'CONVERT-IMPORT > EDF (multiple sampling frequencies)';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(15);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    sProcess.isSeparator = 0;
    % Option: File to import
    sProcess.options.datafile.Comment = 'You will be asked to select EDF/REC files';
    sProcess.options.datafile.Type    = 'label';
    sProcess.options.datafile.Value   = [];
    % Separator
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
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
    {'*.edf;*.rec' , 'EDF Files (*.edf, *.rec)'}, ...
    'Pick EDF files', ...
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

    % Read header
    hdr = ReadEDFHeader(datafiles{iInput});
    
    % If file is not Unique Sampling Rate (USR), convert it
    if numel(unique(hdr.samples))==1 % Unique Sampling Rate (USR) file
        bst_report('Warning', sProcess, sInputs(iInput), sprintf('File [%s] already of unique sampling rate (USR)',[subjectname,ext]))
        usrfilename = datafiles{iInput};
    else % Convert to Unique Sampling Rate (USR) prior to import
        hdr = Convert2USR(hdr);
        usrfilename = hdr.filename;
    end
    
    % Import in BST database
    % Process: Create link to raw file
    sFiles = bst_process(...
        'CallProcess', 'process_import_data_raw', ...
        [], [], ...
        'subjectname', subjectname, ...
        'datafile', {usrfilename, 'EEG-EDF'}, ...
        'channelalign', sProcess.options.channelalign.Value);
    
    OutputFiles{iInput} = sFiles.FileName;
end

end

% -------------------------------------------------------------------------
function hdrUSR = Convert2USR(hdrMSR)
% Convert Multiple Sampling Rates (MSR) EDF file to Unique Sampling Rate
% (USR) EDF file that can be used by brainstorm. 
% 
% "hdrMSR" is the header of MSR file as obtained by <ReadEDFHeader> method
% of <process_sc_import_edf> process. "hdrUSR" is the modified header for
% written USR file.

% Modify header MSR->USR
hdrUSR = hdrMSR;
hdrUSR.filename = strrep(hdrMSR.filename,'.edf','.usr.edf');
% Force integer second duration (problems if, eg, 0.5)
hdrUSR.duration = 1; 
fsMSRMax = max(hdrMSR.samples/hdrMSR.duration); % New unique samprate
hdrUSR.samples(:) = round(fsMSRMax*hdrUSR.duration);
% Adapt number of records
hdrUSR.recordnb = floor(hdrMSR.duration*hdrMSR.recordnb/hdrUSR.duration);
% BST has problems with annotations.. rename them so it wont see it
iAnnot = find(strcmpi(hdrUSR.channelname,'EDF Annotations'));
for i=1:numel(iAnnot)
    hdrUSR.channelname{iAnnot(i)} = sprintf('ANNOTATION%d',i);
end

% Write USR header
[fidUSR,hdrUSR] = process_sc_export_edf('WriteEDFHeader', hdrUSR);

% 
% % % % % minDt = 1/fsMSRMax;
% % % % % TimeBounds = (0:hdrUSR.recordnb-1)*hdrUSR.duration;
% % % % % TimeBounds = [TimeBounds; TimeBounds+hdrUSR.duration-minDt];
TimeBounds = (0:hdrUSR.recordnb-1)*hdrUSR.duration;
TimeBounds = [TimeBounds; TimeBounds+hdrUSR.duration];

try bst_progress('start', 'Converting EDF (MSR->USR)', 'Writing data record...', 0, hdrUSR.recordnb); end
for iRec=1:hdrUSR.recordnb
    try bst_progress('set', iRec); end
    F = ReadEDFData( hdrMSR, TimeBounds(:,iRec));
    if isempty(F), error('End of file overpassed'); end
    process_sc_export_edf('WriteEDFRecord',fidUSR,hdrUSR,F);
end
fclose(fidUSR);

end

% -------------------------------------------------------------------------
function [ hdr ] = ReadEDFHeader( filename )
% (http://www.edfplus.info/specs/edf.html)
% HDR:
% +-- version:        Version number (0)
% +-- patient:        Patient informations
% +-- recording:      Recording informations
% +-- startdate:      Recording date (dd.mm.yy)
% +-- starttime:      Recording time (hh.mm.ss) 
% +-- headersize:     256*(Nc+1)
% +-- reserved:       Reserved field
% +-- recordnb:     * Number of data records
% +-- duration:     * Duration of one data record (sec)
% +-- channelnb:      Number of channels = Nc
% +-- channelname:  * Ncx1 channel labels (e.g. EEG Fpz-Cz or Body temp)
% +-- transducer:     Ncx1 transducer type (e.g. AgAgCl electrode) 
% +-- physdim:      * Ncx1 physical dimension (e.g. uV or degreeC) 
% +-- physmin:      * Ncx1 physical minimum (e.g. -500 or 34) 
% +-- physmax:      * Ncx1 physical maximum (e.g. 500 or 40) 
% +-- digimin:      * Ncx1 digital minimum (e.g. -2048) 
% +-- digimax:      * Ncx1 digital maximum (e.g. 2047) 
% +-- prefilt:        Ncx1 prefiltering (e.g. HP:0.1Hz LP:75Hz)
% +-- samples:      * Ncx1 number of samples in each data record 
% +-- chnreserved:    Ncx1 Reserved field at the end of header
% 
% Author: Jonathan Godbout, 2009-2013

% Open file for reading
fid = fopen(filename,'r','ieee-le');
if fid == -1, error('[read_edf] Cannot open %s.',filename); end
[pathstr, name, format] = fileparts(filename);
if ~ismember(lower(format),{'.edf' '.rec'})
	error('Format %s is not .edf or .rec.',format);
end

hdr.filename = filename;

% Read Header
hdrstr = char(fread(fid,256,'uchar')');
hdr.version = hdrstr(1:8); % ..................................... 8 ASCII
hdr.patient = hdrstr(9:88); % ................................... 80 ASCII
hdr.recording = hdrstr(89:168); % ............................... 80 ASCII
hdr.startdate = hdrstr(169:176); % ............................... 8 ASCII
hdr.starttime = hdrstr(177:184); % ............................... 8 ASCII
hdr.headersize = str2num(hdrstr(185:192)); % ..................... 8 ASCII
hdr.reserved = hdrstr(193:236); % ............................... 44 ASCII
hdr.recordnb = str2num(hdrstr(237:244)); % ....................... 8 ASCII
hdr.duration = str2num(hdrstr(245:252)); % ....................... 8 ASCII
hdr.channelnb = str2num(hdrstr(253:256)); % ...................... 4 ASCII
hdr.channelname = char(fread(fid,[16,hdr.channelnb],'char')');
hdr.transducer = char(fread(fid,[80,hdr.channelnb],'char')');
hdr.physdim = char(fread(fid,[8,hdr.channelnb],'char')');
hdr.physmin = str2num(char(fread(fid,[8,hdr.channelnb],'char')'));
hdr.physmax = str2num(char(fread(fid,[8,hdr.channelnb],'char')'));
hdr.digimin = str2num(char(fread(fid,[8,hdr.channelnb],'char')'));
hdr.digimax = str2num(char(fread(fid,[8,hdr.channelnb],'char')'));
hdr.prefilt = char(fread(fid,[80,hdr.channelnb],'char')');
hdr.samples = str2num(char(fread(fid,[8,hdr.channelnb],'char')'));
hdr.chnreserved = char(fread(fid,[32,hdr.channelnb],'char')');

fclose(fid);

% Compute number of records if unknown
if hdr.recordnb == -1
    f_inf = dir(filename);
    n_bytes_of_data = f_inf.bytes - hdr.length;
    n_int16_of_data = n_bytes_of_data/2;
    hdr.recordnb = n_int16_of_data/sum(hdr.samples);
end

% Convert CHAR matrices to CELL arrays
hdr.channelname = strarray2cell(hdr.channelname);
hdr.transducer  = strarray2cell(hdr.transducer);
hdr.physdim     = strarray2cell(hdr.physdim);
hdr.prefilt     = strarray2cell(hdr.prefilt);
hdr.chnreserved = strarray2cell(hdr.chnreserved);

% EDF type (EDF ou EDF+)
if ismember(upper(hdr.reserved(1:5)),['EDF+C' 'EDF+D'])
    disp('Format EDF+'); %edf_plus = 1; 
    if strcmp(hdr.reserved(1:5), 'EDF+D')
        disp('Data is discontinuous (not supported). Considering as continuous');
    end
end

end
    
function c = strarray2cell(s)

[Nr,Nc] = size(s);
c = cell(Nr,1);
for ii = 1:Nr
    % Remove filling spaces prior
    idx = find(~(s(ii,:)==' '),1,'last');
    if isempty(idx), idx = Nc; end
    c{ii} = s(ii,1:idx);
end

end

% -------------------------------------------------------------------------
function F = ReadEDFData( hdr, TimeBounds, iChannels)
% 
%   F = ReadEDFData( hdr, TimeBounds)             % All Channels
%   F = ReadEDFData( hdr, TimeBounds, iChannels)  % Only channel indexes
% 
% For a single use, "hdr" must have a field "filename" containing the fu

% Open file for reading
fid = fopen(hdr.filename,'r','ieee-le');
if fid == -1, error('[read_edf] Cannot open %s.',hdr.filename); end
[pathstr, name, format] = fileparts(hdr.filename);
if ~ismember(lower(format),{'.edf' '.rec'})
	error('Format %s is not .edf or .rec.',format);
end

if numel(TimeBounds)~=2, error('TimeBounds must be a 2-element vector'); end

if nargin<3, iChannels = 1:hdr.channelnb; end

fs = hdr.samples/hdr.duration;

% Records logistic calculations
rec_length = sum(hdr.samples); % Total number of samples in one record
rec_low  = floor(TimeBounds(1)/hdr.duration); % First record to read
rec_high = ceil(TimeBounds(end)/hdr.duration); % Last record to read
nb_rec_to_read = rec_high - rec_low; % Number of records to read

% Read records
seekf = hdr.headersize + (rec_low*rec_length)*2; % INT16 (2*8bits = 2Bytes)
fseek(fid,seekf,'bof');
try
    data_1d      = fread(fid, rec_length*nb_rec_to_read, 'int16');
    data_records = reshape(data_1d, [rec_length, nb_rec_to_read])';
catch ME
    warning('This end of file was overpassed while reading. No data acquired.')
    F = [];
    return;
end
fclose(fid);

% Extract channels' signal from record form
F = cell(hdr.channelnb,1); % MSR data (cell array for different sized sig)
for ii=1:hdr.channelnb
    % 
    seek = sum(hdr.samples(1:ii-1));
    temp = data_records(:,seek+1:seek+hdr.samples(ii));
    F{ii} = reshape(temp', [1 numel(temp)]);
           
    % Seulement garder la fenêtre d'analyse
    iStart  = round( mod(TimeBounds(1),hdr.duration)*fs(ii) ) + 1;
    iLength = round( (TimeBounds(2)-TimeBounds(1)) * fs(ii) );
    
    F{ii} = F{ii}(iStart:iStart+iLength-1);
end
F = F(iChannels);                   % Take selected channels
F = interpolateLowerSamprate(F);    % Interpolate to have unique samprate
F = digi_to_phys_range(F,hdr,iChannels); % Convert digital values to physical

end

function FUSR = interpolateLowerSamprate(FMSR)

N = cellfun(@numel,FMSR); % Number of samples for each channel
Nmx = max(N);             % Highest samples number
toInterp = find(N<Nmx);   % Interpolate to Nmx all channels < Nmx
for i=1:numel(toInterp)
    FMSR{toInterp(i)} = interpft(FMSR{toInterp(i)},Nmx);
end
FUSR = cell2mat(FMSR); % Convert cell FMSR to matrix FUSR

end

function F = digi_to_phys_range(F,hdr,iChannels)

Nt = size(F,2);

% Change range from analog (uV) to digital values
range_in  = [hdr.digimin(:) hdr.digimax(:)];  % [- +]
range_out = [hdr.physmin(:) hdr.physmax(:)]; % [- +]
range_in  = range_in(iChannels(:),:);
range_out = range_out(iChannels(:),:);

mean_in   = (range_in (:,2) + range_in (:,1))/2 * ones(1,Nt);
mean_out  = (range_out(:,2) + range_out(:,1))/2 * ones(1,Nt);
amp_in    = (range_in (:,2) - range_in (:,1))/2 * ones(1,Nt);
amp_out   = (range_out(:,2) - range_out(:,1))/2 * ones(1,Nt);

F  = (F - mean_in) .* amp_out./amp_in + mean_out;
F  = double(F);

end

