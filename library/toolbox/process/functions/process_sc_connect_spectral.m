function varargout = process_sc_connect_spectral( varargin )
% process_sc_connect_spectral: Computes connectivities from TIMEFREQ data.
% You must compute a time-frequency representation prior to this process.
% 
% For coherency based metrics (Coherence, Magnitude Squared Coherence,
% Imagainary coherence, ...), the common time-frequency representation
% is a spectrogram (Short-Time Fourier Transform).
% 
% For Phase Lag Index (PLI) ans Phase Coherence (or PLV), the common
% time-frequency representation is the Hilbert Transform.
% 
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'CONNECTIVITY > FREQUENCY > Spectral metrics';
    sProcess.FileTag     = 'connect';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1145);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 1;
%     %
%     sProcess.options.N.Comment = 'Number of cross-FFT to mean: ';
%     sProcess.options.N.Type    = 'value';
%     sProcess.options.N.Value   = {10,'',0};
    %
    sProcess.options.ccoh.Comment = 'Complex Coherency | ccoh';
    sProcess.options.ccoh.Type    = 'checkbox';
    sProcess.options.ccoh.Value   = 0;
    %
    sProcess.options.msc.Comment = 'Magnitude Squared Coherence | msc';
    sProcess.options.msc.Type    = 'checkbox';
    sProcess.options.msc.Value   = 0;
    %
    sProcess.options.acoh.Comment = 'Coherence (module of Coherency) | acoh';
    sProcess.options.acoh.Type    = 'checkbox';
    sProcess.options.acoh.Value   = 0;
    %
    sProcess.options.icoh.Comment = 'Imaginary Coherency (Nolte 2004) | icoh';
    sProcess.options.icoh.Type    = 'checkbox';
    sProcess.options.icoh.Value   = 0;
    %
    sProcess.options.ycoh.Comment = 'Absolute Imaginary Coherency | ycoh';
    sProcess.options.ycoh.Type    = 'checkbox';
    sProcess.options.ycoh.Value   = 0;
    %
    sProcess.options.pcoh.Comment = 'Phase Coherence (PLV equivalent) | pcoh';
    sProcess.options.pcoh.Type    = 'checkbox';
    sProcess.options.pcoh.Value   = 0;
%     %
%     sProcess.options.rcoh.Comment = 'rcoh = real(COHERENCY)';
%     sProcess.options.rcoh.Type    = 'checkbox';
%     sProcess.options.rcoh.Value   = 0;
    %
    sProcess.options.upli.Comment = 'Phase Lag Index (Stam 2007) | upli';
    sProcess.options.upli.Type    = 'checkbox';
    sProcess.options.upli.Value   = 0;
%     %
%     sProcess.options.dpli.Comment = 'dpli (directed PLI, to be implemented, Stam 2012)';
%     sProcess.options.dpli.Type    = 'checkbox'; 
%     sProcess.options.dpli.Value   = 0;
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

% Metrics
metrics = [];
if sProcess.options.acoh.Value>0, metrics = [metrics,{'acoh'}]; end
if sProcess.options.icoh.Value>0, metrics = [metrics,{'icoh'}]; end
if sProcess.options.upli.Value>0, metrics = [metrics,{'upli'}]; end
if sProcess.options.ccoh.Value>0, metrics = [metrics,{'ccoh'}]; end
if sProcess.options.msc.Value>0,  metrics = [metrics,{'msc'}]; end
if sProcess.options.ycoh.Value>0, metrics = [metrics,{'ycoh'}]; end
if sProcess.options.pcoh.Value>0, metrics = [metrics,{'pcoh'}]; end
% if sProcess.options.rcoh.Value, metrics = [metrics,{'rcoh'}]; end
% if sProcess.options.dpli.Value, metrics = [metrics,{'dpli'}]; end
Nm = numel(metrics);

Ni = numel(sInputs);
OutputFiles = cell(Ni,Nm);
bst_progress('start', 'Processing', 'Computing Connectivity...', 0, numel(sInputs));
for iInput=1:numel(sInputs)
    bst_progress('set', iInput);
    
    TimefreqMatIn = in_bst_timefreq(sInputs(iInput).FileName, 0);
    
    if size(TimefreqMatIn.TF,1)==numel(TimefreqMatIn.RowNames)
        TimefreqMatIn.TF = Cross(TimefreqMatIn.TF,'product');
    end
    
    if any(ismember(metrics,{'ccoh','msc','acoh','icoh','rcoh','ycoh'}))
        COH = Coherency(TimefreqMatIn.TF);
    end
    if any(ismember(metrics,{'dpli','upli'}))
        PLI = PhaseLagIndex(TimefreqMatIn.TF);
    end
    if any(ismember(metrics,{'pcoh'}))
        PCOH = PhaseCoherence(TimefreqMatIn.TF);
    end
 
    for k = 1:Nm
        filetag = metrics{k};
        TimefreqMat = TimefreqMatIn;
        switch metrics{k}
            case 'ccoh'
                TimefreqMat.TF = COH;        
                TimefreqMat.Method = 'cohere'; % See FIGURE_CONNECT>UpdateColormap
            case 'acoh'
                TimefreqMat.TF = abs(COH);        
                TimefreqMat.Method = 'cohere'; % See FIGURE_CONNECT>UpdateColormap
            case 'msc'
                TimefreqMat.TF = abs(COH).^2;        
                TimefreqMat.Method = 'cohere'; % See FIGURE_CONNECT>UpdateColormap
            case 'rcoh'
                TimefreqMat.TF = real(COH);        
                TimefreqMat.Method = 'corr'; % See FIGURE_CONNECT>UpdateColormap
            case 'icoh'
                TimefreqMat.TF = imag(COH);        
                TimefreqMat.Method = 'corr'; % See FIGURE_CONNECT>UpdateColormap
            case 'ycoh'
                TimefreqMat.TF = abs(imag(COH));        
                TimefreqMat.Method = 'cohere'; % See FIGURE_CONNECT>UpdateColormap
            case 'pcoh'
                TimefreqMat.TF = abs(PCOH);        
                TimefreqMat.Method = 'cohere'; % See FIGURE_CONNECT>UpdateColormap
            case 'upli'
                TimefreqMat.TF = abs(PLI);        
                TimefreqMat.Method = 'corr'; % See FIGURE_CONNECT>UpdateColormap
            case 'dpli'
                TimefreqMat.TF = PLI;        
                TimefreqMat.Method = 'cohere'; % See FIGURE_CONNECT>UpdateColormap
            otherwise
                error('Bad metric')
        end
        TimefreqMat.Comment = [sInputs(iInput).Comment, ' | ', filetag];
%         TimefreqMat.Comment = strrep(sInputs(iInput).Comment,'| spectra' , ['| ', filetag]);
        TimefreqMat.DataType = 'data';
        TimefreqMat.DataFile = sInputs(iInput).DataFile;
        TimefreqMat.Measure = 'other';
        TimefreqMat.Time = [TimefreqMat.Time(1), TimefreqMat.Time(end)];
        TimefreqMat.RowNames = TimefreqMat.RowNames;
        TimefreqMat.RefRowNames = TimefreqMat.RowNames;
        TimefreqMat.Time = [TimefreqMat.Time(1), TimefreqMat.Time(end)];

        TimefreqMat.Options.isMirror = 1;
        TimefreqMat.Options.isSymmetric = 1;
        
        % Group (or not) in frequency bands
        if sProcess.options.isfreqbands.Value   
            if iscell(TimefreqMat.Freqs)
                strMsg = 'Already grouped in frequency bands.';
                bst_report('Warning',sProcess, sInputs(iInput), strMsg);
            else
                % Call function to group by frequency bands
                FreqBands = sProcess.options.freqbands.Value;
                [TimefreqMat] = process_tf_bands('Compute',TimefreqMat, FreqBands, []);
            end
        end

        % Output filename: add file tag
        tagstr = ['timefreq_connectn_',filetag];
        if ~isempty(sInputs(iInput).DataFile)
%             OutputFile = strrep(file_fullpath(sInputs(iInput).DataFile), '.mat', ['_' filetag '.mat']);
%             OutputFile = strrep(OutputFile, [filesep,'data_'], [filesep,'timefreq_connectn_']);
%             OutputFile = file_unique(OutputFile);
            OutputFile = file_fullpath(sInputs(iInput).DataFile);
            OutputFile = strrep(OutputFile, [filesep,'data_'], [filesep,tagstr,'_']);
            OutputFile = file_unique(OutputFile);
        else
            OutputFile = bst_process('GetNewFilename', bst_fileparts(sInputs(iInput).FileName), tagstr);
        end
        % Save file
        bst_save(OutputFile, TimefreqMat, 'v6');
        % Add file to database structure
        db_add_data(sInputs(iInput).iStudy, OutputFile, TimefreqMat);
        OutputFiles{iInput,k} = OutputFile;
    end

end
OutputFiles = OutputFiles(:);

end

function C = Coherency(F)
% F: Nc(Nc+1)/2 x Nt x Nf Cross-Timefreq

% N: number of cross-FFTs to mean for cross-spectra estimation

% % % % % Convolution mean (cross-spectra estimation) with mirror padding
% % % % [Nc,Nt,Nf] = size(F);
% % % % for i=1:Nc
% % % %     for k=1:Nf
% % % %        Fi = Mirror(F(i,:,k),N);
% % % %        Fi = conv(Fi,ones(1,N),'same')/N;
% % % %        F(i,:,k) = Fi(N+1:end-N);
% % % %     end
% % % % end

F = mean(F,2);

[Nx,Nt,Nf] = size(F);
Nc = (-1+sqrt(1+ 8*Nx))/2;

% Normalize cross-spectrum by product of auto-spectra
C = zeros(Nx,Nt,Nf);
for k=1:size(C,2)
    Fk = Vec2Mat(F(:,k,:),'conjugate'); % Convert vector to matrix form
    for i=1:size(Fk,3) % For all frequencies
        d = diag(Fk(:,:,i)); % Auto-terms
        Fk(:,:,i) = Fk(:,:,i)./sqrt(d*d'); % Normalize by products
    end
    Fk = Mat2Vec(Fk); % Reshape Nxx1xNf
    C(:,k,:) = Fk;
end


end

function C = PhaseLagIndex(F)
% F: Nc(Nc+1)/2 x Nt x Nf Cross-Timefreq

p = angle(F);                      % Phase differences
C = mean(sign(imag(exp(1i*p))),2); % PLI (COORAY 2011)

end

function C = PhaseCoherence(F)
% F: Nc(Nc+1)/2 x Nt x Nf Cross-Timefreq

p = angle(F);         % Phase differences
p = exp(1j*p);        % Project to unit circle
C = mean(p,2);        % Circular mean

end

function y = Cross(x,operation)
% Input x is a Nc x Nt x Nf signal (optional 3rd dimension)
% Output is the Nc(Nc+1)/2 x Nt x Nf cross-terms 
% 
% Product:  x(i,:,:) .* conj(x(j,:,:))
% The resulting amplitude is the cross-power and resulting phase is the
% phase difference between each channel combination
% 
% Quotient: x(i,:,:) ./ x(j,:,:)
% The resulting amplitude is the power-ratio and resulting phase is the
% phase difference between each channel combination
% 

if nargin==1, operation = 'product'; end

[Nc,Nt,Nf] = size(x);

Ncx = Nc*(Nc+1)/2;

% Have to adapt to BST 
matidx = find(triu(ones(Nc))>0);
[subi,subj] = ind2sub([Nc,Nc],matidx);

y = zeros(Ncx,Nt,Nf);
switch lower(operation)
    case 'product'
        for i=1:Ncx
            y(i,:,:) = x(subi(i),:,:) .* conj(x(subj(i),:,:));
        end
    case 'quotient'
        for i=1:Ncx
            y(i,:,:) = x(subi(i),:,:) ./ x(subj(i),:,:);
        end
end

% idx = 0;
% y = zeros(Ncx,Nt,Nf);
% switch lower(operation)
%     case 'product'
%         for i=1:Nc
%             for j=i:Nc
%                 idx = idx+1;
%                 y(idx,:,:) = x(i,:,:) .* conj(x(j,:,:));
%             end
%         end
%     case 'quotient'
%         for i=1:Nc
%             for j=i:Nc
%                 idx = idx+1;
%                 y(idx,:,:) = x(i,:,:) ./ (x(j,:,:));
%             end
%         end
% end

end

function TFm = Vec2Mat(TFv,sym)

TFm = sc_bst_connect_format_vec2mat(TFv,sym);

end

function TFv = Mat2Vec(TFm)

TFv = sc_bst_connect_format_mat2vec(TFm);

end


