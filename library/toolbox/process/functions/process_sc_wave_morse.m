function varargout = process_sc_wave_morse( varargin )
% process_sc_wave_morse: Computes Morse continuous analytical wavelet
% transform.
% 
%       OutputFiles = process_sc_wave_morse('Run', sProcess, sInputs)
%       [TF] = process_sc_wave_morse('Morse', x, f, order, moment)
%
% Authors: Jonathan Godbout, 2013 

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'FREQUENCY > Morse Continuous Wavelet Transform';
    sProcess.FileTag     = 'morse';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1050);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    %
    sProcess.options.band.Comment = 'Frequency band: ';
    sProcess.options.band.Type    = 'range';
    sProcess.options.band.Value   = {[11,15],'Hz',2};
    %
    sProcess.options.resolution.Comment = 'Frequency resolution';
    sProcess.options.resolution.Type    = 'value';
    sProcess.options.resolution.Value   = {0.1,'Hz',4};
    %
    sProcess.options.order.Comment = 'Wavelet order';
    sProcess.options.order.Type    = 'value';
    sProcess.options.order.Value   = {20,'',0};
    %
    sProcess.options.moment.Comment = 'Wavelet number of vanishing moments';
    sProcess.options.moment.Type    = 'value';
    sProcess.options.moment.Value   = {10,'',0};
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
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

% Get user parameters
band        = sProcess.options.band.Value{1};
resolution  = sProcess.options.resolution.Value{1};
order       = sProcess.options.order.Value{1};
moment      = sProcess.options.moment.Value{1};

% Loop over the files, compute and save
Ni = numel(sInputs);
OutputFiles = cell(Ni,1);
bst_progress('start', 'Processing', 'Computing Morse...', 0, numel(sInputs));
for iInput=1:numel(sInputs)
    bst_progress('set', iInput);
    
    % Get data and channel structs
    DataMat = in_bst_data(sInputs(iInput).FileName,sInputs(iInput).ChannelFile);
    ChanMat = in_bst_channel(sInputs(iInput).ChannelFile);
        
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
    frequency = band(1):resolution:band(2);
    TF = Morse(x, frequency/fs, order, moment);
     
    Comment = FormatComment(sProcess);

    % SAVE RESULTS --------------------------------------------------------
    TimefreqMat = db_template('timefreqmat');
    TimefreqMat.TF = TF;
    TimefreqMat.Comment = Comment; % Filetag only in file name
%     TimefreqMat.Comment = [FormatComment(sProcess), ' | ', sProcess.FileTag];
    TimefreqMat.DataType = 'data';
    TimefreqMat.DataFile = sInputs(iInput).FileName;
    TimefreqMat.Events = DataMat.Events;
    TimefreqMat.Time = DataMat.Time;
    TimefreqMat.Freqs = frequency;
    TimefreqMat.RowNames = fileChanNames(iChannels);
    TimefreqMat.RefRowNames = [];
    TimefreqMat.Measure = 'none';
    TimefreqMat.Method = sProcess.FileTag;
    TimefreqMat.Options.morse.band = band;
    TimefreqMat.Options.morse.resolution = resolution;
    TimefreqMat.Options.morse.moment = moment;
    TimefreqMat.Options.morse.order = order;
    TimefreqMat.Options.morse.channelNames = channelNames;

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
function Wab = Morse(x,fc,order,moment)
% X is a NchanxNtime data matrix. FC is a vector of Nfreq normalized 
% frequencies (frequencies divided by sampling frequency) at which to 
% compute Morse wavelets coefficients. ORDER and MOMENT are wavelets
% parameters. Returns NchanxNtimexNfreq complex wavelets coefficients.

if nargin<4, moment = 10; end
if nargin<3, order = 20; end

% Add mirror pads to avoid edge effects
Npad = round(10*(1/min(fc))); % Pad: include 10 cycles of lowest frequency
Npad = min(Npad, size(x,2));
x = [x(:,Npad:-1:1), x, x(:,end:-1:end-Npad+1)];

% Sizes of padded signals
[Nc Ns] = size(x);

m = order;
n = moment;

f0 = (((n+0.5)/m).^(1/m))/2/pi; % frequence centrale de l'ondelette de Cauchy

a = f0./fc;
Na = numel(a);

% --- espace de Fourier
wk = 2.*pi/Ns/1*[1:fix(Ns/2)];
wk = [0., wk, -wk(fix((Ns-1)/2):-1:1)];

% --- fft du signal
f = fft(x,[],2);

Wab = zeros(Nc,Ns,Na);
Wab = Wab + 1i*Wab;

for ii=1:Na
%     psi_wk = repmat(Morse_w(a(ii)*wk,n,m).*sqrt(a(ii)),Nc,1);
    psi_wk = repmat(Morse_w(a(ii)*wk,n,m).*a(ii),Nc,1);
    Wab(:,:,ii) = ifft(f.*psi_wk,[],2);
end

% Remove pads
Wab = Wab(:,Npad+1:end-Npad,:);


end

function y = Morse_w(w,n,m)
y  = zeros(size(w));
wp = w(w>0);
y(w>0) = wp.^n.*exp(-wp.^m);
end

