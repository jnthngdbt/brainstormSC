function varargout = process_sc_export_edf( varargin )
% process_sc_export_edf: Export data to EDF
%  
%       hdr = process_sc_export_edf('WriteEDFHeader',hdr)
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'EXPORT > EDF';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1745);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    % Records duration
    sProcess.options.duration.Comment = 'Records (epochs) duration: ';
    sProcess.options.duration.Type    = 'value';
    sProcess.options.duration.Value   = {1, 'sec', 0};
    % Mode
    sProcess.options.mode.Comment = {'Write new EDF file (header + data)','Write only a header file (.edfhdr)'};
    sProcess.options.mode.Type    = 'radio';
    sProcess.options.mode.Value   = 1;
    % Downsample
    sProcess.options.downsample.Comment = 'Downsampling factor: ';
    sProcess.options.downsample.Type    = 'value';
    sProcess.options.downsample.Value   = {1, '', 0};
    % === 
    sProcess.options.label1.Comment = [ ...
'Notes: ' ...
'<BR> --- sampling frequency is divided by this value' ...
'<BR> --- must be a factor (divisor) of the sampling frequency (otherwise, data may corrupt)' ...
'<BR> --- implicit low-pass filtering at 75% of new Nyquist frequency to avoid aliasing' ...
        ];
    sProcess.options.label1.Type    = 'label';
    % Options: Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor name(s) or type(s) (Empty=All): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = [];
    % === 
    sProcess.options.label2.Comment = [ ...
        '<HTML>If you specify types instead of names (or left empty): ' ...
'<BR> --- channel list may not be the same for files having different channel files' ...
'<BR>  If you only specify names:' ...
'<BR> --- files having unavailable channel(s) will be skipped' ...
        ];
    sProcess.options.label2.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

savepath = [uigetdir, filesep]; % Export path
Ni = numel(sInputs);

% Parameters
downSample = sProcess.options.downsample.Value{1};
if downSample<1, downSample = 1; end

% Write data
for iInput=1:Ni

    % Get data structs
    sData = in_bst_data(sInputs(iInput).FileName);
    ChanMat = in_bst_channel(sInputs(iInput).ChannelFile);
    
    fsIn = 1/(sData.Time(2)-sData.Time(1)); % Sampling frequency
    fsOut = fsIn/downSample; % Sampling frequency after downsampling
    
    % Channel names in current file
    fileChanNames = {ChanMat.Channel.Name};
    % Get channels
    [iChannels, channelNames, WrnMsg, ErrMsg] = sc_bst_channel_find(ChanMat.Channel, sProcess.options.sensortypes.Value);
    % Report if necessary
    if ~isempty(WrnMsg)
        bst_report('Warning', sProcess, sInputs(iInput), WrnMsg);
    end
    if ~isempty(ErrMsg)
        bst_report('Error', sProcess, sInputs(iInput), ErrMsg);
        continue; % Skip this file
    end

    % === HEADER ===
    hdr.channelname = fileChanNames(iChannels);
    hdr.channelnb   = numel(iChannels);
    hdr.startdate   = '01.01.01';
    hdr.starttime   = '00.00.00';
    hdr.duration    = sProcess.options.duration.Value{1};
    hdr.recordnb    = floor(sData.Time(end)/hdr.duration);
    hdr.gain        = 10^6*ones(hdr.channelnb,1); % BST works in V, we want uV
    hdr.physmin     = -2^14*ones(hdr.channelnb,1); % 16384, gives 0.5uV resolution in EDF
    hdr.physmax     = 2^14*ones(hdr.channelnb,1)-1;
    hdr.physdim     = cell(1,hdr.channelnb);
    hdr.physdim(:)  = {'uV'}; % Micro symbol not supported
    hdr.digimin     = -32768*ones(hdr.channelnb,1);
    hdr.digimax     = 32767*ones(hdr.channelnb,1);
    hdr.samples     = round(hdr.duration*fsOut)*ones(hdr.channelnb,1);
        
    % Redneck units checking...
    if mean(hdr.physmax)<10^2 % Assumed to be Volts
        strMsg = sprintf('Signals in [%s] seems to be Volts. Scale to have uV.',savename);
        bst_report('Error', sProcess, sInputs, strMsg);
        continue;
    end 
    
    % Create and write the header
    savename = [savepath,sInputs(iInput).SubjectName,'.EXPORT.edf'];
    switch sProcess.options.mode.Value
        case 2 % Write only a header file (.edfhdr)
            hdr.filename = [savename,'hdr'];
            [fid,hdr] = WriteEDFHeader(hdr);
            if fid==-1 % If open fails, report and pass to next file
                strMsg = sprintf('Could not write file [%s]',hdr.filename);
                bst_report('Error', sProcess, sInputs, strMsg);
                continue;
            end
            fclose(fid);
            continue % Done. Go to next file
        otherwise
            hdr.filename = savename;
            [fid,hdr] = WriteEDFHeader(hdr);
            if fid==-1 % If open fails, report and pass to next file
                strMsg = sprintf('Could not write file [%s]',hdr.filename);
                bst_report('Error', sProcess, sInputs, strMsg);
                continue;
            end
    end
    
    % === WRITE DATA ===
    recLength = round(hdr.duration*fsIn); % Length of record to read (BEFORE downsampling)
% % % % % %     SamplesBounds = (0:hdr.recordnb-1)*recLength + 1; 
    SamplesBounds = (0:hdr.recordnb-1)*recLength; % BST read seems to be 0-based indexing
    SamplesBounds = [SamplesBounds; SamplesBounds + recLength - 1];
    bst_progress('start', sprintf('Exporting EDF (%d/%d)',iInput,Ni), 'Writing data record...', 0, hdr.recordnb);
    for iRec = 1:hdr.recordnb
        bst_progress('set', iRec);
        % Get window data
        switch upper(sData.F.format)
            case {'EEG-EDF', 'EEG-BDF'}
                ChannelRange = [min(iChannels),max(iChannels)];
                [F, TimeVector] = in_fread(sData.F, 1, SamplesBounds(:,iRec),ChannelRange);
                F = F(iChannels-min(iChannels)+1,:);
            case {'CTF', 'CTF-CONTINUOUS'}
                ChannelRange = min(iChannels):max(iChannels);
                [F, TimeVector] = in_fread(sData.F, 1, SamplesBounds(:,iRec),ChannelRange);
                F = F(iChannels-min(iChannels)+1,:);
            otherwise
                [F, TimeVector] = in_fread(sData.F, 1, SamplesBounds(:,iRec), iChannels);
        end
        % Downsampling
        if downSample>1 % Do it
            NyquistFreq = fsOut/2;
            LowPass = (3/4)*NyquistFreq; % 3/4 of Nyquist frequency after downsampling
            F = process_bandpass('Compute', F, fsIn, 0, LowPass);
            F = F(:,1:downSample:end); % DOWNSAMPLE
        end
        
        WriteEDFRecord(fid,hdr,F);
    end
        
    fclose(fid);
end

% % % % % [F, TimeVector] = in_fread(sFile, 1, SamplesBounds, iChannels);

% === EXIT ===
OutputFiles = {sInputs.FileName}; % Return input file names

end

% ------------------------------------------------------------------------
function [fid,hdr] = WriteEDFHeader(hdr)
% http://www.edfplus.info/specs/edf.html
% EVERY NUMBER MUST BE AN INTEGER (e.g.: duration must be an integer number
% of seconds)
% 
% HDR: * = absolutely required
% +-- filename:     * Full path and name of EDF file to write
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

fid = fopen(hdr.filename, 'wb','ieee-le');
if fid==-1, return; end % If open fails, abort

% Fill missing and detect errors
if ~isfield(hdr,'filename'),    error('Missing FILENAME'); end
if ~isfield(hdr,'version'),     hdr.version = 0; end
if ~isfield(hdr,'patient'),     hdr.patient = '<no patient info>'; end
if ~isfield(hdr,'recording'),   hdr.recording = '<no local recording info>'; end
if ~isfield(hdr,'startdate'),   hdr.startdate = '01.01.01'; end
if ~isfield(hdr,'starttime'),   hdr.starttime = '00.00.00'; end
if ~isfield(hdr,'headersize'),  hdr.headersize = 256*(1+numel(hdr.channelname)); end
if ~isfield(hdr,'reserved'),    hdr.reserved = ''; end
if ~isfield(hdr,'recordnb'),    error('Missing RECORDNB'); end
if ~isfield(hdr,'duration'),    error('Missing DURATION'); end
if ~isfield(hdr,'channelnb'),   hdr.channelnb = numel(hdr.channelname); end
if ~isfield(hdr,'channelname'), error('Missing CHANNELNAME'); end
if ~isfield(hdr,'transducer'),  hdr.transducer = cell(hdr.channelnb,1); end
if ~isfield(hdr,'physdim'),     error('Missing PHYSDIM'); end
if ~isfield(hdr,'physmin'),     error('Missing PHYSMIN'); end
if ~isfield(hdr,'physmax'),     error('Missing PHYSMAX'); end
if ~isfield(hdr,'digimin'),     error('Missing DIGIMIN'); end
if ~isfield(hdr,'digimax'),     error('Missing DIGIMAX'); end
if ~isfield(hdr,'prefilt'),     hdr.prefilt = cell(hdr.channelnb,1); end
if ~isfield(hdr,'samples'),     error('Missing SAMPLES'); end
if ~isfield(hdr,'chnreserved'), hdr.chnreserved = cell(hdr.channelnb,1); end

% Convert header fields to fixed size char arrays and WRITE in file
hdrFields = ...
{'version'    ,'patient'   ,'recording','startdate','starttime' , ...
 'headersize' ,'reserved'  ,'recordnb' ,'duration' ,'channelnb' , ...
 'channelname','transducer','physdim'  ,'physmin'  ,'physmax'   , ...
 'digimin'    ,'digimax'   ,'prefilt'  ,'samples'  ,'chnreserved'};
hdrFieldsStrSize = ...
[8            ,80          ,80         ,8          ,8           , ...
 8            ,44          ,8          ,8          ,4           , ...
 16           ,80          ,8          ,8          ,8           , ...
 8            ,8           ,80         ,8          ,32          ];
for i=1:numel(hdrFields)
    hdrstr = convert2fixedstringarray(hdr.(hdrFields{i}), hdrFieldsStrSize(i));
    fwrite(fid, hdrstr', 'uchar');
end

end

function s = convert2fixedstringarray(c,maxlength)

if isempty(c), c = {[]}; end

if ischar(c), Nc = size(c,1);
else          Nc = numel(c);
end
s = []; % String matrix to be
for ii=1:Nc
    if iscell(c)
        L = sprintf(sprintf('%%-%ds',maxlength),c{ii});
    elseif isnumeric(c)
        L = sprintf(sprintf('%%-%dd',maxlength),c(ii));
    elseif ischar(c)
        L = sprintf(sprintf('%%-%ds',maxlength),c(ii,:));
    end
    end_i = min(numel(L),maxlength);
    s = [s; L(1:end_i)];
end

end

% ------------------------------------------------------------------------
function WriteEDFRecord(fid,hdr,F)

Nt = size(F,2);

% Apply gain to convert to uV
if isfield(hdr,'gain')
    F = (hdr.gain*ones(1,Nt)) .* F;
end

% Change range from analog (uV) to digital values
range_in  = [hdr.physmin(:) hdr.physmax(:)]; % [- +]
range_out = [hdr.digimin(:) hdr.digimax(:)];  % [- +]

mean_in  = (range_in (:,2) + range_in (:,1))/2 * ones(1,Nt);
mean_out = (range_out(:,2) + range_out(:,1))/2 * ones(1,Nt);

amp_in  = (range_in (:,2) - range_in (:,1))/2 * ones(1,Nt);
amp_out = (range_out(:,2) - range_out(:,1))/2 * ones(1,Nt);

F = (F - mean_in) .* amp_out./amp_in + mean_out;

fwrite(fid,int16(F)','int16'); % VECTORIZED writing (by records)

end

