function varargout = process_sc_connect_partialcorr( varargin )
% process_sc_connect_partialcorr: Computes partial correlation from data.
% Based on Marrelec 2006.
% 
% External call
% 
%   C = process_sc_connect_partialcorr('PartialCorrelation',x)
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'CONNECTIVITY > Partial Correlation';
    sProcess.FileTag     = 'pcorr';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1135);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
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
    % === IS FREQ BANDS
    sProcess.options.sep2.Type    = 'separator';
    sProcess.options.isfreqbands.Comment = 'Compute in frequency bands (name/freqs/function):';
    sProcess.options.isfreqbands.Type    = 'checkbox';
    sProcess.options.isfreqbands.Value   = 1;
    % === FREQ BANDS
    sProcess.options.freqbands.Comment = '';
    sProcess.options.freqbands.Type    = 'groupbands';
    sProcess.options.freqbands.Value   = bst_get('DefaultFreqBands');
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

Ni = numel(sInputs);
OutputFiles = cell(Ni,1);
bst_progress('start', 'Processing', 'Computing Partial Correlation...', 0, numel(sInputs));
for iInput=1:numel(sInputs)
    bst_progress('set', iInput);
    
    % Get data and channel structs
    DataMat = in_bst_data(sInputs(iInput).FileName,sInputs(iInput).ChannelFile);
    ChanMat = in_bst_channel(sInputs(iInput).ChannelFile);
       
    % Sampling frequency
    fs = 1/(DataMat.Time(2)-DataMat.Time(1));
    
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

    nChannels = numel(iChannels);
    % Compute (or not) in frequency bands
    if sProcess.options.isfreqbands.Value   
        Freqs = sProcess.options.freqbands.Value;
        FreqBounds = process_tf_bands('GetBounds',Freqs);
        nFreqs = size(FreqBounds,1);
        C = zeros(nChannels*(nChannels+1)/2, 1, nFreqs);
        for iFreq = 1:nFreqs
            x = bst_bandpass(DataMat.F(iChannels,:),fs,FreqBounds(iFreq,1),FreqBounds(iFreq,2));
            C(:,:,iFreq) = PartialCorrelation(x);
        end
    else
        x = DataMat.F(iChannels,:);
        C = PartialCorrelation(x);
        Freqs = 0;
    end

    % Create output timefreq connectivity structure
    filetag = sProcess.FileTag;
    TimefreqMat = db_template('timefreqmat');
    TimefreqMat.TF = C;        
    TimefreqMat.Method = 'corr'; % See FIGURE_CONNECT>UpdateColormap
    TimefreqMat.Comment = [sInputs(iInput).Comment, ' | ', filetag];
    TimefreqMat.DataType = 'data';
    TimefreqMat.DataFile = sInputs(iInput).FileName;
    TimefreqMat.Measure = 'other';
    TimefreqMat.Time = [DataMat.Time(1), DataMat.Time(end)];
    TimefreqMat.Freqs = Freqs; 
    TimefreqMat.RowNames = fileChanNames(iChannels);
    TimefreqMat.RefRowNames = TimefreqMat.RowNames;
    TimefreqMat.Options.isMirror = 1;
    TimefreqMat.Options.isSymmetric = 1;

    % Save in database
    OutputFile = file_fullpath(sInputs(iInput).FileName);
    OutputFile = strrep(OutputFile, [filesep,'data_'], [filesep,'timefreq_connectn_',filetag]);
    OutputFile = file_unique(OutputFile);
%     if ~isempty(sInputs(iInput).DataFile)
%         OutputFile = file_fullpath(sInputs(iInput).DataFile);
%         OutputFile = strrep(OutputFile, [filesep,'data_'], [filesep,tagstr,'_']);
%         OutputFile = file_unique(OutputFile);
%     else
%         OutputFile = bst_process('GetNewFilename', bst_fileparts(sInputs(iInput).FileName), tagstr);
%     end
    % Save file
    bst_save(OutputFile, TimefreqMat, 'v6');
    % Add file to database structure
    db_add_data(sInputs(iInput).iStudy, OutputFile, TimefreqMat);
    OutputFiles{iInput} = OutputFile;
    
end

end

function C = PartialCorrelation(x)
% Based on MARRELEC 2006
% Input X is a Nc x Nt x Nf (3rd dimension is optional)
% Output C is [Nc(Nc+1)/2 x 1 x Nf] connectivity

[Nc,Nt,Nf] = size(x);
C = zeros(Nc,Nc,Nf);
for i=1:Nf
    K = cov(x(:,:,i)'); % Covariance matrix
    Y = pinv(K);        % Inversion via pseudo-inverse
    % Normalize cross-term with auto-terms product
    d = diag(Y);
    Ci = -Y./sqrt(d*d');
    Ci(eye(Nc)==1) = 1; % Set diagonal to 1 (by definition)
    C(:,:,i) = Ci;
end
% Convert in vector form
C = sc_bst_connect_format_mat2vec(C);

end

