function varargout = process_sc_spectrogram( varargin )
% PROCESS_SC_SPECTROGRAM: Computes spectrogram (serie of FFTs), also known
% as Short-Time Fourier Transform.
% 
%       OutputFiles = process_sc_spectrogram('Run', sProcess, sInputs)
%       [F,freq,t] = process_sc_spectrogram('Spectrogram', x, cfg)
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'FREQUENCY > Spectrogram (Short-Time Fourier Transform)';
    sProcess.FileTag     = 'spectra';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1045);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    %
    sProcess.options.duration.Comment = 'FFT Length: ';
    sProcess.options.duration.Type    = 'value';
    sProcess.options.duration.Value   = {4,'sec',2};
    %
    sProcess.options.band.Comment = 'Frequency limits: ';
    sProcess.options.band.Type    = 'range';
    sProcess.options.band.Value   = {[0,50],'Hz',2};
    %
    sProcess.options.resolution.Comment = 'Frequency resolution (0-> 1/[FFT Length]): ';
    sProcess.options.resolution.Type    = 'value';
    sProcess.options.resolution.Value   = {[],'Hz',2};
    %
    sProcess.options.overlap.Comment = 'Overlap: ';
    sProcess.options.overlap.Type    = 'value';
    sProcess.options.overlap.Value   = {50,'%',0};
    %
    sProcess.options.window.Comment = 'Tapering window: ';
    sProcess.options.window.Type    = 'combobox';
    sProcess.options.window.Value   = {4,{'rectangular','cosine','hamming','hann'}};
    %
    sProcess.options.preprocess.Comment = 'Preprocessing: ';
    sProcess.options.preprocess.Type    = 'combobox';
    sProcess.options.preprocess.Value   = {1,{'none','demean','detrend'}};
    % Options: Sensor types
    sProcess.options.sensortypes.Comment = 'Sensor name(s) or type(s) (Empty=All): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = [];
    % === 
    sProcess.options.label2.Comment = [ ...
        '<HTML>If you specify types instead of names (or left empty): ' ...
'<BR> --- channel list may not be the same for files having different channel files' ...
'<BR> --- this may be problematic if consistency across subjects (or conditions) is needed' ...
'<BR> If you only specify names:' ...
'<BR> --- files having unavailable channel(s) will be skipped' ...
        ];
    sProcess.options.label2.Type    = 'label';
    %     % === OVERWRITE
% %     sProcess.options.overwrite.Comment = 'Overwrite input files';
% %     sProcess.options.overwrite.Type    = 'checkbox';
% %     sProcess.options.overwrite.Value   = 1;                                             
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
%     Comment = sProcess.Comment;
    b = sProcess.options.band.Value{1};
    o = sProcess.options.overlap.Value{1};
    d = sProcess.options.duration.Value{1};
    r = sProcess.options.resolution.Value{1};
    if ~isempty(r) && r>0
    else
        r = 1/d;
    end
    Comment = sprintf('Spectrogram [%.2f-%.2fHz] [%.2fHz] [%.2fs] [%d%%]',b(1),b(2),r,d,o);
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

% Get user parameters
band        = sProcess.options.band.Value{1};
overlap     = sProcess.options.overlap.Value{1};
duration    = sProcess.options.duration.Value{1};
window      = sProcess.options.window.Value{2}{sProcess.options.window.Value{1}};
preprocess  = sProcess.options.preprocess.Value{2}{sProcess.options.preprocess.Value{1}};
resolution  = sProcess.options.resolution.Value{1};

if ~isempty(resolution) && resolution>0
else
    resolution = 1/duration;
end

% Loop over the files, compute and save
Ni = numel(sInputs);
OutputFiles = cell(Ni,1);
bst_progress('start', 'Processing', 'Computing Spectrogram...', 0, numel(sInputs));
for iInput=1:numel(sInputs)
    bst_progress('set', iInput);
    
    % Get data and channel structs
    DataMat = in_bst_data(sInputs(iInput).FileName,sInputs(iInput).ChannelFile);
    ChanMat = in_bst_channel(sInputs(iInput).ChannelFile);
       
    % Make sure that data is long enough for at least one FFT
    if (DataMat.Time(end)-DataMat.Time(1)) < duration
        msgstr = 'Data is too short';
        bst_report('Error', sProcess, sInputs(iInput), msgstr);
        continue;
    end
    
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
    
    % COMPUTE -------------------------------------------------------------
    x = DataMat.F(iChannels,:);
    fs = 1/(DataMat.Time(2)-DataMat.Time(1)); 	% Sampling frequency
    cfg.band        = band/fs;                  % Frequency band
    cfg.length      = round(duration*fs);       % Segments length
    cfg.overlap     = overlap/100;              % Segments overlap
    cfg.window      = window;                   % Tapering window name
    cfg.preprocess  = preprocess;               % Segments data preprocess
    cfg.nfft        = round((1/resolution)*fs); % FFT computation length
    [TF,F,T] = Spectrogram(x, cfg);
    T = T/fs + DataMat.Time(1);
    F = F*fs;
     
    Comment = FormatComment(sProcess);
    if size(TF,2)==1
        T = [DataMat.Time(1),DataMat.Time(end)];
        Comment = strrep(Comment,sprintf('[%d%%]',overlap), '[0%]');
%         Comment = strrep(Comment,sprintf('[%d%%]',overlap), '[--%]');
    end

    % SAVE RESULTS --------------------------------------------------------
    TimefreqMat = db_template('timefreqmat');
    TimefreqMat.TF = TF;
    TimefreqMat.Comment = Comment; % Filetag only in file name
%     TimefreqMat.Comment = [FormatComment(sProcess), ' | ', sProcess.FileTag];
    TimefreqMat.DataType = 'data';
    TimefreqMat.DataFile = sInputs(iInput).FileName;
    TimefreqMat.Events = DataMat.Events;
    TimefreqMat.Time = T;
    TimefreqMat.Freqs = F;
    TimefreqMat.RowNames = fileChanNames(iChannels);
    TimefreqMat.RefRowNames = [];
    TimefreqMat.Measure = 'none';
    TimefreqMat.Method = sProcess.FileTag;
    TimefreqMat.Options.spectra.band = band;
    TimefreqMat.Options.spectra.overlap = overlap;
    TimefreqMat.Options.spectra.duration = duration;
    TimefreqMat.Options.spectra.window = window;
    TimefreqMat.Options.spectra.preprocess = preprocess;
    TimefreqMat.Options.spectra.resolution = resolution;
    TimefreqMat.Options.spectra.channelNames = channelNames;

    % Output filename: add file tag
    OutputFile = strrep(file_fullpath(sInputs(iInput).FileName), '.mat', ['_' sProcess.FileTag '.mat']);
    OutputFile = strrep(OutputFile, [filesep,'data_'], [filesep,'timefreq_']);
    OutputFile = file_unique(OutputFile);
    % Save file
    bst_save(OutputFile, TimefreqMat, 'v6');
    % Add file to database structure
    db_add_data(sInputs(iInput).iStudy, OutputFile, TimefreqMat);
    OutputFiles{iInput} = OutputFile;

end


end

% ========================================================================
function [F,freq,t] = Spectrogram(x,cfg)
% Computes a spectrogram (or Short-Time Fourier Transform) on a NcxNt data
% matrix 'x'. 'cfg' is a configuration structure. The default values for 
% 'cfg' (or if it is not sent as input) corresponds to a standard
% normalized FFT. Valid fields are:
% 
% cfg.length:   number of data samples for each single FFT
% cfg.window: 	name of the tapering window ('rectangular', 'hann',
%               'hamming' or 'cosine')
% cfg.overlap:  overlap between FFTs. If <1, considered as a fraction of
%               FFT length. If >1, considered as number of samples.
% cfg.preprocess: name of pre-processing ('none', 'demean' or 'detrend')
% cfg.nfft:     number of points on which to compute FFT
% cfg.band:     limiting frequency band (fraction of sampling frequency)

[Nc,Nt] = size(x);          % Number of time samples Nt and channels Nc

% Set defaults (if nothing specified, defaults equals a normalized FFT)
if nargin<2, cfg = struct; end
if ~isfield(cfg,'length'),      cfg.length = Nt; end
if ~isfield(cfg,'window'),      cfg.window = 'rectangular'; end
if ~isfield(cfg,'overlap'),     cfg.overlap = 0; end
if ~isfield(cfg,'preprocess'),  cfg.preprocess = 'none'; end
if ~isfield(cfg,'nfft'),        cfg.nfft = cfg.length; end
if ~isfield(cfg,'band'),        cfg.band = []; end

% If 0<=overlap<1, convert into number of samples (fraction of length)
if cfg.overlap<1, cfg.overlap = round(cfg.overlap*cfg.length); end

Nl = cfg.length;                % Number of samples for FFT
Nf = floor(cfg.nfft/2);      	% Number of frequencies
No = cfg.overlap;               % Number of overlaping samples 
Ns = Nl-No;                     % Number of samples steps between windows
Nw = floor((Nt-Nl)/Ns)+1;       % Number of windows (FFTs)

% TAPERING WINDOW
switch lower(cfg.window)
    case 'hann',        w = 0.50-0.50*cos(2*pi*(0:Nl-1)/(Nl-1));
    case 'hamming',     w = 0.54-0.46*cos(2*pi*(0:Nl-1)/(Nl-1));
    case 'cosine',      w = sin(pi*(0:Nl-1)/(Nl-1));
    case 'rectangular', w = ones(1,Nl);
    otherwise,          error('Invalid tapering window name');
end
w = ones(Nc,1)*w;       % Convert to matrix (duplicate to each channel)

% PREPROCESSING
switch lower(cfg.preprocess)
    case 'none',        hP = @(x) x;                        % None
    case 'demean',     	hP = @(x) detrend(x','constant')';    % De-mean
    case 'detrend',    	hP = @(x) detrend(x','linear')';      % De-trend
    otherwise,          error('Invalid preprocessing name');
end

F = complex(zeros(Nc,Nw,Nf));
t = zeros(1,Nw); % Time vector (samples) (center of windows)
for iWin=1:Nw
    start_i = (iWin-1)*Ns+1;            % Window start sample
    end_i = start_i+Nl-1;               % Window end sample
    t(iWin) = start_i+floor(Nl/2)-1;    % Window center position
    xi = x(:, start_i:end_i);           % Window signal
    xi = hP(xi);                        % Pre-processing (apply handle)
    xi = w .* xi;                       % Tapered data window
    Xi = fft(xi',cfg.nfft)'/Nl;       	% FFT of data window
    F(:,iWin,:) = Xi(:,1:Nf);           % Stack FFT in output matrix
end

% Frequency vector (0Hz to Nyquist) (fraction of sampling frequency)
freq = (0:cfg.nfft)/cfg.nfft;
freq = freq(1:Nf);

% If a frequency band is specficied, return only this band
if ~isempty(cfg.band) && numel(cfg.band)==2
    band_i = (freq >= cfg.band(1) & freq <= cfg.band(2));
    freq = freq(band_i);
    F = F(:,:,band_i);
end

end

