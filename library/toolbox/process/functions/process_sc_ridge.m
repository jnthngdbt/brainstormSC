function varargout = process_sc_ridge( varargin )
% process_sc_ridge: Extract ridge (oscillator) from a time-frequency
% representation using local maxima (frequency dimension) in specified
% bands.
% 
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'FREQUENCY > Ridge (oscillator extraction)';
    sProcess.FileTag     = 'ridge';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1070);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
    % === IS FREQ BANDS
    sProcess.options.sep2.Type    = 'separator';
    sProcess.options.isfreqbands.Comment = 'Group by frequency bands (name/freqs/function):';
    sProcess.options.isfreqbands.Type    = 'checkbox';
    sProcess.options.isfreqbands.Value   = 1;
    % === FREQ BANDS
    sProcess.options.freqbands.Comment = '';
    sProcess.options.freqbands.Type    = 'groupbands';
    sProcess.options.freqbands.Value   = bst_get('DefaultFreqBands');
end

%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
%     N = sProcess.options.N.Value{1};
%     Comment = sprintf('Complex Coherency [%d]',N);
%     Comment = sprintf('Coherency');
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>



filetag = sProcess.FileTag;
Ni = numel(sInputs);
OutputFiles = cell(Ni,1);
bst_progress('start', 'Processing', 'Extracting ridge...', 0, numel(sInputs));
for iInput=1:numel(sInputs)
    bst_progress('set', iInput);
    
    TimefreqMatIn = in_bst_timefreq(sInputs(iInput).FileName, 0);

    % Frequency bands
    if sProcess.options.isfreqbands.Value    
        FreqBands = process_tf_bands('GetBounds',sProcess.options.freqbands.Value);
    else
        FreqBands = [TimefreqMatIn.Freqs(1), TimefreqMatIn.Freqs(end)];
    end
    
    [TFr, iFreq, ridgeMap] = Ridge(TimefreqMatIn.TF,TimefreqMatIn.Freqs,FreqBands);
    
    TimefreqMat = TimefreqMatIn;
    TimefreqMat.TF = TFr;
    TimefreqMat.Freqs = sProcess.options.freqbands.Value;
    TimefreqMat.Comment = [sInputs(iInput).Comment, ' | ', filetag];
    TimefreqMat.DataType = 'data';
    TimefreqMat.DataFile = sInputs(iInput).DataFile;
    TimefreqMat.Measure = 'none';
    

    % Output filename: add file tag
    tagstr = ['timefreq_',filetag];

    OutputFile = bst_process('GetNewFilename', bst_fileparts(sInputs(iInput).FileName), tagstr);
    % Save file
    bst_save(OutputFile, TimefreqMat, 'v6');
    % Add file to database structure
    db_add_data(sInputs(iInput).iStudy, OutputFile, TimefreqMat);
    OutputFiles{iInput} = OutputFile;
end
OutputFiles = OutputFiles(:);

end

% =========================================================================
function [TFr, iFreq, ridgeMap] = Ridge(TF,Freq,FreqBands)
% TF is a Nchan x Ntime x Nfreq timefrequency matrix. FREQ is the vector of
% the Nfreq frequency domain values. FreqBands is a Nband x 2 matrix of
% limiting sub-bands in which to separately compute the ridge. 
% 
% Output TFR is the Nchan x Ntime x Nband ridge matrix. IFREQ is a Nchan x
% Ntime x Nband frequency indexes of the frequency maxima (indexes of
% FREQ).

if ~exist('FreqBands','var'), FreqBands = [Freq(1), Freq(end)]; end

[Nc,Nt,Nf] = size(TF);
Nb = size(FreqBands,1); % Number of bands

TFr = zeros(Nc,Nt,Nb);
iFreq = zeros(Nc,Nt,Nb);
RidgeMap = cell(1,Nb);
for iBand=1:Nb
    iFreqBand = find(Freq>=FreqBands(iBand,1) & Freq<=FreqBands(iBand,2));
    if ~isempty(iFreqBand)
        [TFr(:,:,iBand), iFreq(:,:,iBand), ridgeMap{iBand}] = ...
            RidgeBand(TF(:,:,iFreqBand));
        % Convert back subband freq indexes to full freq domain indexes
        iFreq(:,:,iBand) = iFreqBand(iFreq(:,:,iBand));
    end
end

end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [TFr, iFreq, ridgeMap] = RidgeBand(TF)

[Nc,Nt,Nf] = size(TF);
                    
TFa = abs(TF);
[mmx,iFreq] = max(TFa,[],3);

TFr = zeros(Nc,Nt);
for cc=1:Nc
    for ss=1:Nt
        TFr(cc,ss) = TF(cc, ss, iFreq(cc,ss));
    end
end

ridgeMap = [];
if nargout>=3
    i1 = ones(Nt,1)*(1:Nc); i2 = (1:Nt)'*ones(1,Nc); i3 = iFreq';
    ridgeMap = zeros(Nc,Nt,Nf); ridgeMap(sub2ind(size(ridgeMap),i1(:),i2(:),i3(:))) = 1;
    ridgeMap = squeeze(sum(ridgeMap))';
end

end

