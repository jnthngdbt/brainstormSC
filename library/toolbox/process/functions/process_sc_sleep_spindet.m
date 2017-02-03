function varargout = process_sc_sleep_spindet( varargin )
% process_sc_sleep_spindet: Detect spindles in raw files. Returns raw files
% with new events.
%
% Author: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>

%     evtList = {'REM','NREM1','WAKE'};

    % Description the process
    sProcess.Comment     = 'SLEEP > Spindles Detection';
    sProcess.FileTag     = '| spindet';
    sProcess.Description = '<HTML>Detect sleep spindles in RAW data';
    sProcess.Category    = 'custom';
    sProcess.Index       = sc_bst_process_index(2045);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    %
    sProcess.options.label0.Comment = '<HTML>(*) = Important parameters. (e) For expert supervised classification';
    sProcess.options.label0.Type    = 'label';
    % ---------------------------------- ~1min / 20min of recording, with 10 channels (Macbook Pro 2011)
    sProcess.options.label1.Comment = '<HTML>========== ANALYSIS STEP 1: Events creation (adaptative segmentation using wavelet ridges)';
    sProcess.options.label1.Type    = 'label';
    % 
    sProcess.options.channel.Comment = '(*) Detection channels name/type (coherent spindle activity): ';
    sProcess.options.channel.Type    = 'text';
    sProcess.options.channel.Value   = 'f3,f4,c3,c4,p3,p4';
    % 
    sProcess.options.channeltopo.Comment = 'Additionnal channels to include for topography computation: ';
    sProcess.options.channeltopo.Type    = 'text';
    sProcess.options.channeltopo.Value   = '';
    % 
    sProcess.options.freqBand.Comment = '(*) Spindle frequency band: ';
    sProcess.options.freqBand.Type    = 'range';
    sProcess.options.freqBand.Value   = {[10,16],'Hz',1};
    %
    sProcess.options.freqRes.Comment = 'Frequency resolution: ';
    sProcess.options.freqRes.Type    = 'value';
    sProcess.options.freqRes.Value   = {0.1,'Hz',3};
%     %
%     sProcess.options.moment.Comment = 'Morse wavelet vanishing moments: ';
%     sProcess.options.moment.Type    = 'value';
%     sProcess.options.moment.Value   = {10,'',0};
%     %
%     sProcess.options.order.Comment = 'Morse wavelet order: ';
%     sProcess.options.order.Type    = 'value';
%     sProcess.options.order.Value   = {20,'',0};
    % ----------------------------------
    sProcess.options.label2.Comment = '<HTML>========== ANALYSIS STEP 2: Candidates selection (threshold applied on amplitude and sigma-index features)';
    sProcess.options.label2.Type    = 'label';
    % 
    sProcess.options.evtrejectin.Comment = 'Events created inside those markers are artifacts (empty=None): ';
    sProcess.options.evtrejectin.Type    = 'text';
    sProcess.options.evtrejectin.Value   = 'harmact, artefact, muscle';
    % 
    sProcess.options.evtrejectout.Comment = 'Events created outside those markers are artifacts (empty=None): ';
    sProcess.options.evtrejectout.Type    = 'text';
    sProcess.options.evtrejectout.Value   = 'Sample Section';
    % 
    sProcess.options.evth0.Comment = '(*) Events created inside those markers are used to construct Null hypothesis ("no spindle" model) (empty=All): ';
    sProcess.options.evth0.Type    = 'text';
    sProcess.options.evth0.Value   = 'REM';
    %
    sProcess.options.pvalue.Comment = '(*) Threshold (p-value rejecting Null hypothesis): ';
    sProcess.options.pvalue.Type    = 'value';
    sProcess.options.pvalue.Value   = {0.06,'',3};
%     %
%     sProcess.options.interval.Comment = 'Merging interval: ';
%     sProcess.options.interval.Type    = 'value';
%     sProcess.options.interval.Value   = {1,'ms',0};
%     %
%     sProcess.options.duration.Comment = 'Time window around events for topography computation: ';
%     sProcess.options.duration.Type    = 'value';
%     sProcess.options.duration.Value   = {0.5,'ms',0};
%     %
%     sProcess.options.topoFreqBand.Comment = 'Frequency band for topography computation: ';
%     sProcess.options.topoFreqBand.Type    = 'range';
%     sProcess.options.topoFreqBand.Value   = {[10,16],'Hz',1};
    % ----------------------------------
    sProcess.options.label3.Comment = '<HTML>========== ANALYSIS STEP 3: Dendrogram (hierarchy of distances on frequency and spatial features)';
    sProcess.options.label3.Type    = 'label';
    % 
    sProcess.options.evtexclout.Comment = 'Only consider candidates inside those markers (empty=All): ';
    sProcess.options.evtexclout.Type    = 'text';
    sProcess.options.evtexclout.Value   = 'Sample Section, NREM2, NREM3, NREM4';
%     %
%     sProcess.options.metric.Comment = 'Hierarchical clustering metric (see PDIST): ';
%     sProcess.options.metric.Type    = 'combobox';
%     sProcess.options.metric.Value   = {1,{'euclidean'}};
    %
    sProcess.options.method.Comment = 'Hierarchical clustering linkage method (see LINKAGE): ';
    sProcess.options.method.Type    = 'combobox';
    sProcess.options.method.Value   = {1,{'average','centroid','complete','median','single','ward','weighted'}};
    % ----------------------------------
    sProcess.options.label4.Comment = '<HTML>========== ANALYSIS STEP 4: Classification (supervised if expert events available)';
    sProcess.options.label4.Type    = 'label';
%     %
%     sProcess.options.criterion.Comment = 'Classification cut-off criterion (see CLUSTER): ';
%     sProcess.options.criterion.Type    = 'combobox';
%     sProcess.options.criterion.Value   = {4,{'uniform','maxclust','inconsistent','distance'}};
    %
    sProcess.options.cutoff.Comment = 'Dendrogram distance cut-off (classification): ';
    sProcess.options.cutoff.Type    = 'value';
    sProcess.options.cutoff.Value   = {0.25,'',3}; % for 11-16Hz range, normalized 0.25 distance is ~1.25Hz
    % 
    sProcess.options.evtexpert.Comment = '(e) Those markers are expert spindles (for supervised class selection) (empty=None): ';
    sProcess.options.evtexpert.Type    = 'text';
    sProcess.options.evtexpert.Value   = '';
    % 
% % %     sProcess.options.expertcutoff.Comment = '(e) Minimum ratio of expert events to recover from class selection: ';
    sProcess.options.expertcutoff.Comment = '(e) Minimum ratio of events concordant with expert events to recover from class selection: ';
    sProcess.options.expertcutoff.Type    = 'value';
    sProcess.options.expertcutoff.Value   = {80,'%',0};
    % ----------------------------------
    sProcess.options.separator.Type = 'separator';
    sProcess.options.separator.Comment = ' ';
    %
    sProcess.options.tasks.Comment = {'Start at step 1','Start at step 2','Start at step 3','Start at step 4','Only process results'};
    sProcess.options.tasks.Type    = 'radio';
    sProcess.options.tasks.Value   = 1;
    % === 
    sProcess.options.isfigres.Comment = 'Show results figures';
    sProcess.options.isfigres.Type    = 'checkbox';
    sProcess.options.isfigres.Value   = 1;
    % === 
    sProcess.options.isviewsds.Comment = 'Show analysis parameters structure';
    sProcess.options.isviewsds.Type    = 'checkbox';
    sProcess.options.isviewsds.Value   = 0;
    % === 
    sProcess.options.isexportcsv.Comment = 'Export markers in CSV file';
    sProcess.options.isexportcsv.Type    = 'checkbox';
    sProcess.options.isexportcsv.Value   = 0;

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

    sProtocol = bst_get('ProtocolInfo');

    bst_progress('start', 'Spindle Detection', 'Initiating');
    for iFile = 1:length(sInputs)

        % ===== GET DATA =====
        % Load the descriptors
        DataMat = in_bst_data(sInputs(iFile).FileName);
        sChan = in_bst_channel(sInputs(iFile).ChannelFile);
        sFile = DataMat.F;
        
        % Channel business 
        iChanTopo = []; % All channels
        if ~isempty(sProcess.options.channeltopo.Value)
            iChanTopo = channel_find(sChan.Channel, sProcess.options.channeltopo.Value); % Cautious with CHANNEL_FIND; it returns a sorted list (may not be the order of how it is specified)
        end
        iChanSpin = channel_find(sChan.Channel, sProcess.options.channel.Value); % Cautious with CHANNEL_FIND; it returns a sorted list (may not be the order of how it is specified)
        iChanTopo = union(iChanTopo,iChanSpin);
        % Remove bad channels
        iBadChan = find(DataMat.ChannelFlag~=1);
        iChanSpin = setdiff(iChanSpin,iBadChan);
        iChanTopo = setdiff(iChanTopo,iBadChan);
                
        % ===== SPINDLE DETECTION SCRIPT (SDS) =====
        
% % % %         SDS_FILENAME = strrep(file_fullpath(sInputs(iFile).FileName),'.mat','_sds.mat');
% % % %         SDS_FILENAME = strrep(SDS_FILENAME,'0raw_','');
        SDS_FILENAME = [sProtocol.STUDIES,filesep,sInputs(iFile).SubjectName,'_sds.mat'];

        if exist(SDS_FILENAME,'file')>0, 
            load(SDS_FILENAME); 
        else
            sProcess.options.tasks.Value = 1; % Force run entire script
        end % sds
        
        % Add or update content of raw file (contains channels, events...)
        sds.file = sFile;

        % =================================================================
        % ======= ANALYSIS STEP 1
        % =================================================================
        if sProcess.options.tasks.Value == 1
            % PARAMETERS
            sds.channel.idx = iChanTopo;
            sds.channel.name = {sChan.Channel(iChanTopo).Name};
            sds.channel.loc = [sChan.Channel(iChanTopo).Loc];
            sds.channel.spindle.idx = iChanSpin;
            sds.channel.spindle.name = {sChan.Channel(iChanSpin).Name};
            sds.channel.spindle.loc =  [sChan.Channel(iChanSpin).Loc];
            sds.frequency.resolution = sProcess.options.freqRes.Value{1};
            sds.frequency.sigma = sProcess.options.freqBand.Value{1};
            sds.frequency.low = [4,10]; % Huupponen 2007 [4,10]
            sds.frequency.high = [20,40]; % Huupponen 2007 [20,40]
            sds.frequency.domain  = sds.frequency.low(1):sds.frequency.resolution:sds.frequency.high(2);
            sds.frequency.lowIdx = find(sds.frequency.domain >= sds.frequency.low(1) & sds.frequency.domain <= sds.frequency.low(2));
            sds.frequency.highIdx = find(sds.frequency.domain >= sds.frequency.high(1) & sds.frequency.domain <= sds.frequency.high(2));
            sds.frequency.sigmaIdx = find(sds.frequency.domain >= sds.frequency.sigma(1) & sds.frequency.domain <= sds.frequency.sigma(2));
            sds.wavelet.moment = 10; % sProcess.options.moment.Value{1};
            sds.wavelet.order  = 20; % sProcess.options.order.Value{1};
            % CHECKS
            if any(sum(sds.channel.loc)==0)
                strMsg = 'Some selected channels have invalid locations. Check channel file.';
                bst_report('Error', sProcess, sInputs(iFile), strMsg);
                continue
            end
            % METHODS
            bst_progress('start', 'Spindle Detection', 'Setting Projector');
            sds = ProjectorSet(sds);
            bst_progress('start', 'Spindle Detection', 'Adaptative segmentation');
            sds = Segment(sds); 
            % FINISH
            save(SDS_FILENAME,'sds');
        end

        % =================================================================
        % ======= ANALYSIS STEP 2
        % =================================================================
        % At this stage, adaptative segmentation of entire file is done. We
        % have ~70000 markers for an 8 hour night. Here we do the following
        % FLAG MARKERS + THRESHOLD CANDIDATES
        if sProcess.options.tasks.Value <= 2
            % PARAMETERS
            sds.flag.h0_l        = ExtractEvtNames(sProcess.options.evth0.Value);
            sds.flag.rejectin_l  = ExtractEvtNames(sProcess.options.evtrejectin.Value);
            sds.flag.rejectout_l = ExtractEvtNames(sProcess.options.evtrejectout.Value);
            sds.threshold.pvalue   = sProcess.options.pvalue.Value{1};
            sds.threshold.feature     = {'amplitude','sigmaindex'};
            sds.threshold.onlyhighforsigma = 0; 
% % %             sds.merge.interval    	= sProcess.options.interval.Value{1}; % Seconds
% % %             sds.merge.feature      	= 'amplitude'; % Feature to consider for retained event
% % %             sds.topography.duration	= sProcess.options.duration.Value{1};
            % METHODS
            bst_progress('start', 'Spindle Detection', 'Thresholding');
%             sds = FlagMarkers(sds); if ischar(sds), bst_report('Error',sProcess,sInputs(iFile),sds); continue; end
            sds = FlagMarkersH0(sds); if ischar(sds), bst_report('Error',sProcess,sInputs(iFile),sds); continue; end
            sds = Threshold(sds);     if ischar(sds), bst_report('Error',sProcess,sInputs(iFile),sds); continue; end
            % FINISH
            save(SDS_FILENAME,'sds');
        end

        % =================================================================
        % ======= ANALYSIS STEP 3
        % =================================================================
        % At this stage, candidates markers passing the null hypothesis
        % p-value threshold are identified and all features are extracted.
        % Here we do the following:
        % AGGLOMERATIVE HIERARCHICAL CLUSTERING
        if sProcess.options.tasks.Value <= 3
            % PARAMETERS
            sds.flag.excludeout_l = ExtractEvtNames(sProcess.options.evtexclout.Value);
            sds.hierarchy.feature = {'frequency'        , 'median'};
            sds.hierarchy.bin     = {sds.frequency.domain(sds.frequency.sigmaIdx), -0.1:0.05:0.1};
            sds.hierarchy.metric  = 'euclidean'; %sProcess.options.metric.Value{2}{sProcess.options.metric.Value{1}};
            sds.hierarchy.method  = sProcess.options.method.Value{2}{sProcess.options.method.Value{1}};
            % METHODS
            bst_progress('start', 'Spindle Detection', 'Hierarchical clustering');
            sds = FlagMarkersCandidate(sds);
            sds = Hierarchy(sds);
            % FINISH
            save(SDS_FILENAME,'sds');
        end
        
        % =================================================================
        % ======= ANALYSIS STEP 4
        % =================================================================
        % CLASSIFICATION + IDENTIFICATION
        if sProcess.options.tasks.Value <= 4
            % PARAMETERS
            sds.flag.expert_l    = ExtractEvtNames(sProcess.options.evtexpert.Value);
            sds.class.criterion    = 'distance'; %sProcess.options.criterion.Value{2}{sProcess.options.criterion.Value{1}};
            sds.class.cutoff       = sProcess.options.cutoff.Value{1};
            sds.class.expertcutoff = sProcess.options.expertcutoff.Value{1};
            % METHODS
            bst_progress('start', 'Spindle Detection', 'Classification');
            sds = Classification(sds);
            bst_progress('start', 'Spindle Detection', 'Merging near events of same class');
            sds = Merge(sds);
            bst_progress('start', 'Spindle Detection', 'Representing and identifying classes');
            sds = Representation(sds);
            sds = Identification(sds);
            % FINISH
            save(SDS_FILENAME,'sds');
            
        end
        
        % =================================================================
        % ======= RESULTS PROCESSING
        % =================================================================
        
        % Show main figure with multiple subplots
        if sProcess.options.isfigres.Value>0
            ShowMainFigure(sds,sInputs(iFile).FileName)
        end
        
        % Show analysis parameters structure
        if sProcess.options.isviewsds.Value
            view_struct(rmfield(sds,'file'));
        end
        
        % Export a CSV file of spindle markers
        if sProcess.options.isexportcsv.Value>0
            bst_progress('start', 'Spindle Detection', 'Exporting events to CSV');
            ExportCSV(sds,strrep(SDS_FILENAME,'.mat','.csv'));
        end
        
        
        bst_progress('start', 'Spindle Detection', 'Creating BST events');
        sFile = ClearEvents(sFile);
        sEvent = CreateEvents(sds);
        
        
        % ===== SAVE RESULT =====
        % Only save changes if something was detected
        if ~isempty(sEvent)
            % Update events structure
            if ~isfield(sFile, 'events') || isempty(sFile.events)
                sFile.events = sEvent;
            else
                sFile.events = [sFile.events, sEvent];
            end
            % Save updated file
            if ~isequal(DataMat.F, sFile)
                DataMat.F = sFile;
                % Keep only time window
                DataMat.Time = DataMat.Time([1, end]);
                % Save file definition
                save(file_fullpath(sInputs(iFile).FileName), '-struct', 'DataMat');
            end
        else
            bst_report('Warning', sProcess, sInputs(iFile), ['No event detected.']);
        end
    end
    % Return all the input files
    OutputFiles = {sInputs(iFile).FileName};
            
    bst_progress('stop');
end

% =========================================================================
% ========== ANALYSIS STEP 1
% =========================================================================

% -------------------------------------------------------------------------
function sds = ProjectorSet(sds)
% Compute a "Spindles" projector from which a mixing matrix is computed
% from first component. By default, a projector is built from an array of
% 1; the resulting mixing matrix corresponds to the mean of selected
% channels. 
% 
% Further developpments could explore the computation of a
% projector from some spindles marked by the expert.

% Build mask vector of selected channels
sds.channel.spindle.mask = zeros(numel(sds.file.channelmat.Channel),1);
sds.channel.spindle.mask(sds.channel.spindle.idx) = 1;

% Build BST projector from artificial data: single sample topography with 1
% at selected channels and 0 elsewhere. There is only 1 component that
% explains 100% of this topo. Seems heavy to do this, but at least the
% structure is set for further developments (adaptative projector from
% expert events)
topoModel = ones(numel(sds.channel.spindle.idx),1);
sds.channel.spindle.projector = process_ssp(...
    'Compute',topoModel,sds.channel.spindle.mask);

% First component (Ui'*Ui = sum(Ui.^2) = 1)
Ui = sds.channel.spindle.projector.Components(sds.channel.spindle.idx,1);
% Mixing matrix
sds.channel.spindle.mix = Ui';

% % % % % % Mixing matrix for 1st component
% % % % % M = Ui*Ui';
% % % % % % 
% % % % % sds.channel.spindle.mix = M(1,:);

end

% -------------------------------------------------------------------------
function sds = Segment(sds)
% Morse ridge adaptative segmentation

fs = sds.file.prop.sfreq;

% We analyse by window blocks for memory (cannot load 8hr nights). To avoid
% segmentation edge effects, windows have overlaps. Segments positioned
% inside edges will be ignored
seg.window.length = 2^14;
seg.window.edge = round(5*fs); 
seg.window.step = seg.window.length - 2*seg.window.edge;

% Windows sample bounds to extract from raw
seg.window.bounds = 1: seg.window.step: sds.file.prop.samples(end);
seg.window.bounds = [seg.window.bounds; seg.window.bounds+seg.window.length-1];
seg.window.bounds(2,end) = sds.file.prop.samples(end);
if diff(seg.window.bounds(:,end)) < 2*seg.window.edge
    seg.window.bounds(:,end) = [];
end

seg.window.epoch = 1;

seg.marker = struct;

% Loop and compute on all window blocks
nWin = size(seg.window.bounds,2);
bst_progress('start', 'Segmenting', 'Windowed segmentation', 0, nWin);
for iWin=1:nWin
    bst_progress('set', iWin);
    
    % ============ TIME_FREQUENCY BUSINESS
    % Get window data
    switch upper(sds.file.format)
        case {'EEG-EDF', 'EEG-BDF'}
            ChannelRange = [min(sds.channel.idx),max(sds.channel.idx)];
            [F, TimeVector] = in_fread(sds.file, ...
                seg.window.epoch, ...
                seg.window.bounds(:,iWin), ...
                ChannelRange);
            F = F(sds.channel.idx-min(sds.channel.idx)+1,:);
        otherwise
            [F, TimeVector] = in_fread(sds.file, ...
                seg.window.epoch, ...
                seg.window.bounds(:,iWin), ...
                sds.channel.idx);
    end
    % Extract spindle channels signals
    iChanSpinRel = cell2mat(arrayfun(@(x) find(x==sds.channel.idx), sds.channel.spindle.idx,'UniformOutput',0));
    Fspn = F(iChanSpinRel,:);
    % Apply projector (mixing matrix)
    Fm = sds.channel.spindle.mix*Fspn;
    % Morse wavelet time-frequency representation (broadband)
    TF = process_sc_wave_morse('Morse', Fm, ...
        sds.frequency.domain/fs, ...
        sds.wavelet.order, ...
        sds.wavelet.moment);
    % Ridge in sigma band
    [TFr, iFreq, ridgeMap] = process_sc_ridge('Ridge', TF, ...
        sds.frequency.domain, ...
        sds.frequency.sigma);
    % ============ SEGMENT MARKERS BUSINESS
    % Get ridge local maxima (temporal dimension) -> segments format
    [mrk] = SegmentMaxima(abs(TFr));
    % Remove segments in window edges
    if iWin~=1,    mrk = sc_class_marker('Remove',mrk, mrk.position_i<seg.window.edge); end
    if iWin~=nWin, mrk = sc_class_marker('Remove',mrk, mrk.position_i>(numel(Fm)-seg.window.edge)); end
    mrk = rmfield(mrk,'channel_i'); % Not valid here
    mrk = rmfield(mrk,'frequency_i'); % Not valid here
    % Add topography feature
    mrk = SegmentTopo(mrk,F,sds);
    % Convert some segments features to context
    mrk.frequency   = sds.frequency.domain(iFreq(mrk.position_i));
    % Locate maxima positions as the max of the real signal
    [mrk] = MaximaRegister(real(TFr),fs,mrk);
    % Add other frequency features
    mrk.low   = mean(abs(TF(:, mrk.position_i, sds.frequency.lowIdx)), 3);
    mrk.high  = mean(abs(TF(:, mrk.position_i, sds.frequency.highIdx)), 3);
    mrk.sigma = mean(abs(TF(:, mrk.position_i, sds.frequency.sigmaIdx)), 3);
    % Convert local positions to full file positions
    mrk.position_i  = mrk.position_i + seg.window.bounds(1,iWin) - 1;

    % Add window segments markers to main markers container
    seg.marker = sc_class_marker('Add',seg.marker, mrk);
    
end

% Rename maxima
seg.marker.amplitude = seg.marker.maxima;
seg.marker = rmfield(seg.marker,'maxima');
% Add sigma index (TODO: ADD CHECKBOXES TO CONSIDER OR NOT LOW AND HIGH BANDS)
seg.marker.sigmaindex = 2*seg.marker.amplitude ./ (seg.marker.low + seg.marker.high);
% Add position (sec)
seg.marker.position = seg.marker.position_i/fs;

% Add to main structure
sds.segment = seg;

end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [mrk] = SegmentMaxima(x)
% Maxima operator. X is a Nchan x Ntime x Nfreq data matrix.

% If X is complex, the local maxima is defined as the maximum of real(X)
% around maxima of abs(X) (TODO)

[Nc,Nt,Nf] = size(x);

mrk = [];
for i=1:Nc
    for j=1:Nf
        mrk = maxima_ij(x(i,:,j),mrk,i,j);
    end
end

end

function [mrk] = maxima_ij(x,mrk,cc,ss)

if isempty(mrk)
    mrk.channel_i = [];
    mrk.position_i = [];
    mrk.maxima = [];
    mrk.frequency_i = [];
end

dx = diff(x);
dx = [dx(1), dx];

% y = zeros(size(x));
position_i = (dx(1:end-1).*dx(2:end))<0; % Find zeros-crossings
position_i = position_i & dx(1:end-1)>0; % Maxima
position_i = find(position_i);
% y(position_i)    = x(position_i);
mrk.position_i      = [mrk.position_i, position_i];
mrk.maxima          = [mrk.maxima, x(position_i)];
mrk.channel_i       = [mrk.channel_i, cc*ones(size(position_i))];
mrk.frequency_i     = [mrk.frequency_i, ss*ones(size(position_i))];

end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [mrk] = SegmentTopo(mrk,F,sds)

FS = sds.file.prop.sfreq;

% Duration (samples) of data window for topo computation
duration_i = round(0.5*FS);
% % % % % duration_i = round(sds.topography.duration*sds.file.prop.sfreq);

% Filtering options
LowPass = sds.frequency.sigma(2);
HighPass = sds.frequency.sigma(1);

% Samples bounds around markers
nTime = size(F,2);
SampleBounds = mrk.position_i - round(duration_i/2);
SampleBounds = [SampleBounds; (mrk.position_i + round(duration_i/2))];
SampleBounds(SampleBounds<1) = 1;
SampleBounds(SampleBounds>nTime) = nTime;

nMrk = size(SampleBounds,2); % Number of markers
nChan = size(F,1); % Number of channels

mrk.topography = zeros(nChan,nMrk);
mrk.coordinates = zeros(3,nMrk);
for iMrk=1:nMrk
    % Extract data window
    Fi = F(:,SampleBounds(1,iMrk):SampleBounds(2,iMrk));
    % Filter in specific band
    Fi = detrend(Fi','linear')';
    Fi = bst_bandpass(Fi, FS, HighPass, LowPass);
    % SVD: topographic components (U vectors)
    [U,S,V] = svd(Fi);
    % Topography is the first component
    mrk.topography(:,iMrk) = U(:,1);
end
% Get location at maximum of topography
% Row 1: front(+)-back(-)
% Row 2: left(+)-right(-) (weird but yes)
% Row 3: up(+)-down(-)
[maxTopo,idxTopo] = max(abs(mrk.topography));
mrk.coordinates = sds.channel.loc(:,idxTopo); 
% Associated channel label
mrk.channel = sds.channel.name(idxTopo);

% For an entire night and a big montage, could be very big
mrk = rmfield(mrk,'topography');

% Median position (front(+)-back(-))
mrk.median = mrk.coordinates(1,:);
% Right-left position (left(+)-right(-)) (weird but yes)
mrk.rightleft = mrk.coordinates(2,:);
end

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
function [mrk] = MaximaRegister(x,fs,mrk)
% MRK must have position_i and frequency features. It will relocate the
% maxima position_i to the local maxima inside a 1/frequency_i window
% around original position_i.
% 
% X must be one signal (1 x Ntime)

N = numel(mrk.position_i);    % Number of markers
Ns = numel(x); % Number of samples
HL = round(0.5*fs./mrk.frequency); % Half-lengths
for i=1:N
    pi = mrk.position_i(i);
    iStart = max(pi-HL(i),1);
    iEnd = min(pi+HL(i),Ns);
    [vMax,iMax] = max(x(iStart:iEnd));
    mrk.position_i(i) = mrk.position_i(i)+iMax-HL(i)-1;
end
mrk.position_i(mrk.position_i<=0) = 1;

end

% =========================================================================
% ========== ANALYSIS STEP 2
% =========================================================================

% -------------------------------------------------------------------------
function [sds] = FlagMarkersH0(sds)

mrk = sds.segment.marker; % All segmentation markers

nMrk = numel(mrk.position_i); % Number of markers

% Flag to reject markers (artefacts) outside some markers
if isempty(sds.flag.rejectout_l)
    mrk.rejectout_l = false(1,nMrk); % Default mode = None
else
    % Create a marker from file event
    ref = FlagMarkersBuild(sds, sds.flag.rejectout_l);
    % Flag segments with that marker
    mrk = sc_class_marker('Flag', mrk, ref, 'second', 'rejectout_l');
    mrk.rejectout_l = ~(mrk.rejectout_l>0);
end

% Flag to reject markers (artefacts) inside some markers
if isempty(sds.flag.rejectin_l)
    mrk.rejectin_l = false(1,nMrk); % Default mode = None
else
    % Create a marker from file event
    ref = FlagMarkersBuild(sds, sds.flag.rejectin_l);
    % Flag segments with that marker
    mrk = sc_class_marker('Flag', mrk, ref, 'second', 'rejectin_l');
    mrk.rejectin_l = mrk.rejectin_l>0;
end

% Flag events at the edge of frequency band
mrk.freqedge_l = ismember(mrk.frequency, sds.frequency.sigma);

% Null hypothesis markers flag
if isempty(sds.flag.h0_l)
    mrk.h0_l = true(1,nMrk); % Default mode = All
else
    % Create a marker from file event
    ref = FlagMarkersBuild(sds, sds.flag.h0_l);
    % Flag segments with that marker
    mrk = sc_class_marker('Flag', mrk, ref, 'second', 'h0_l');
    mrk.h0_l = mrk.h0_l>0;
end

% Check if ok
if ~any(mrk.h0_l)
    errMsg = 'Could not find specified null hypothesis events!';
    sds = errMsg;
    return;
end

% Post-process H0 flags to reject artifacts
% Reject artifacts
mrk.h0_l = mrk.h0_l & ~mrk.rejectin_l;
mrk.h0_l = mrk.h0_l & ~mrk.rejectout_l;
% Do not consider frequency band edge effect
mrk.h0_l = mrk.h0_l & ~mrk.freqedge_l;

% Check if ok
if ~any(mrk.h0_l)
    errMsg = 'All null hypothesis events have been rejected!';
    sds = errMsg;
    return;
end

% OUTPUT
sds.segment.marker = mrk;

end

% -------------------------------------------------------------------------
function sds = Threshold(sds)

% INPUTS
featName = sds.threshold.feature; % Names of features to threshold
pval = sds.threshold.pvalue; % Threshold: p-value on null (H0)
mrk = sds.segment.marker;
Nf = length(featName);
Nm = numel(mrk.position_i);
Nh = sum(mrk.h0_l);

% Recompute on the fly sigma index
if sds.threshold.onlyhighforsigma>0
    mrk.sigmaindex = mrk.amplitude ./ mrk.high;
else
    mrk.sigmaindex = 2*mrk.amplitude ./ (mrk.low + mrk.high);
end

% For each feature, compute threshold: p-value on H0 CDF
featThres = zeros(Nf,1);
mrk.threshold_l = true(1,Nm);
for ii=1:Nf
    % Get H0 values
    featH0Val = mrk.(featName{ii})(mrk.h0_l); % Get H0 values
    % Get value at percentile p 
    featH0Val = sort(featH0Val,'descend');
    featThres(ii) = featH0Val(round(pval*numel(featH0Val)));
    % Old way of getting this value
% % % % % %     [featPDF featPDFBin] = hist(featH0Val,Nh); % PDF
% % % % % %     featCDF = cumsum(featPDF)/Nh;  % CDF
% % % % % %     featThres(ii) = featPDFBin(find(featCDF>=(1-pval),1)); % Alpha of CDF

    % Threshold-pass flag feature
    mrk.threshold_l = mrk.threshold_l & mrk.(featName{ii})>featThres(ii); % Flag
end

sds.threshold.value = featThres;
sds.segment.marker = mrk;

end

% =========================================================================
% ========== ANALYSIS STEP 3
% =========================================================================

% -------------------------------------------------------------------------
function [sds] = FlagMarkersCandidate(sds)

mrk = sds.segment.marker; % All segmentation markers

nMrk = numel(mrk.position_i); % Number of markers

% Flag to exlude markers outside some markers (at the end of detection)
if isempty(sds.flag.excludeout_l)
    mrk.excludeout_l = false(1,nMrk); % Default mode = None
else
    % Create a marker from file event
    ref = FlagMarkersBuild(sds, sds.flag.excludeout_l);
    % Flag segments with that marker
    mrk = sc_class_marker('Flag', mrk, ref, 'second', 'excludeout_l');
    mrk.excludeout_l = ~(mrk.excludeout_l>0);
end

% Candidate flag: respect threshold, no artifact, no exclusions
mrk.candidate_l = mrk.threshold_l;
mrk.candidate_l = mrk.candidate_l & ~mrk.rejectin_l;
mrk.candidate_l = mrk.candidate_l & ~mrk.rejectout_l;
mrk.candidate_l = mrk.candidate_l & ~mrk.excludeout_l;
mrk.candidate_l = mrk.candidate_l & ~mrk.freqedge_l;

% OUTPUT
sds.segment.marker = mrk;
sds.candidate.marker = sc_class_marker('Get',mrk,mrk.candidate_l>0);

end

% -------------------------------------------------------------------------
function sds = Hierarchy(sds)

% Inputs
mrk = sds.candidate.marker; % Candidate markers
metric = sds.hierarchy.metric;
method = sds.hierarchy.method;
featBin = sds.hierarchy.bin;

% Extract matrix of specified features used for classification
featVal = sc_class_marker('Matrix',mrk,sds.hierarchy.feature);
[Nf,Ne] = size(featVal);

% Normalize each feature for them to range in [0;1]
featValMin = zeros(Nf,Ne);
featValMax = zeros(Nf,Ne);
for iFeat=1:Nf
    featValMin(iFeat,:) = min(featBin{iFeat})*ones(1,Ne);
    featValMax(iFeat,:) = max(featBin{iFeat})*ones(1,Ne);
end
featValNorm = (featVal-featValMin)./(featValMax-featValMin);

% Create hierarchy (statistics toolbox)
dist = pdist(featValNorm',metric);  % Pair-wise distances
dist = dist+0.0001*rand(size(dist)); %%%%%%%%% LINKAGE FREEZES IF MANY DIST VALUES ARE EQUAL; ADD SMALL RANDOMNESS.....
link = linkage(dist,method);        % Hierarchy tree
clearvars dist; % Really big
hFig = figure; % To catch the dendrogram that will print
[tmpb,tmpa,rank] = dendrogram(link,Ne);   % Tree dendrogram ordering
close(hFig); % Close the catched dendrogram

% Compute 2D histogram of both features (MUST HAVE 2 FEATURES FOR CLASS)
[featHist,C] = hist3(featVal', featBin);

mrk.classrank = rank(:)';

% Output
sds.candidate.marker = mrk;
sds.hierarchy.vector = featVal;
% % sds.hierarchy.dist = dist; % DONT SAVE IT! Too big
sds.hierarchy.link = link;
sds.hierarchy.rank = rank;
sds.hierarchy.hist = featHist;

end

% =========================================================================
% ========== ANALYSIS STEP 4
% =========================================================================

% -------------------------------------------------------------------------
function sds = Classification(sds)

% Inputs
mrk = sds.candidate.marker; % Candidate markers
criterion = sds.class.criterion;
cutoff = sds.class.cutoff;
link = sds.hierarchy.link;

% Clustering: cut the tree
switch lower(criterion)
    case 'uniform'
        [link,num] = ClusterUniform(link,cutoff); % Home made
    case 'maxclust'
        num = cluster(link,'maxclust',cutoff);
    otherwise
        num = cluster(link,'cutoff',cutoff,'criterion',criterion);
end

% Add classification features to candidate markers
mrk.classnum = num(:)';

% Output
sds.candidate.marker = mrk;
sds.hierarchy.class = mrk.classnum;

end

function [Z T] = ClusterUniform(data, ratio, method)
% Classification using CLUSTER. Number of classes is such that the ratio of
% the second biggest to the biggest class is lower than RATIO.
% 
%   [Z T] = UNIFORM(Z,RATIO)
%   [Z T] = UNIFORM(D,RATIO,METHOD)

if nargin==3
    D = data;
    Z = linkage(D, method);
elseif nargin==2
    Z = data;
end 
    
for ii=2:1000
    disp(['Trying with ' num2str(ii) ' classes..'])
    T = cluster(Z, 'maxclust', ii); % Classify
    % Size of all classes
    N = [];
    Tj = unique(T);
    for jj=1:numel(Tj)
        N = [N sum(T==Tj(jj))];
    end
    % Size ratio of 2 biggest classes
    N = sort(N,'descend');
    r = N(2)/N(1);
    % Return if class size uniformity condition is met 
    if r>=ratio, return; end
end

end

% -------------------------------------------------------------------------
function sds = Merge(sds)
% Merge "near" events of same class.

% INPUTS
mrk  = sds.candidate.marker; % Candidate markers to merge
nMrk = numel(mrk.position);  % Number of markers
% CONSTANT
dt = 0.5;                    % Merging interval (s)

% Loop over classes
nClass = max(mrk.classnum);
mrk.merge_l = true(1,nMrk); % We initialy keep everything
for iClass=1:nClass
    idx = mrk.classnum==iClass;
    mrkClass = sc_class_marker('Get', mrk, idx);
    mrkClass = MergeClass(mrkClass,dt);
    mrk.merge_l(idx) = mrkClass.merge_l;
end

% OUTPUTS
sds.candidate.marker = mrk;
sds.candidate.merge.marker = sc_class_marker('Get', mrk, mrk.merge_l>0);
sds.candidate.merge.interval = dt;

end

function mrk = MergeClass(mrk,dt)
% Merge "near" events of same class.

% INPUTS
feat = mrk.amplitude; % Feature to select merged events representant
pos  = mrk.position;  % Markers positions (s)
nMrk = numel(pos); % Number of markers

% Find those markers that are too near
x = [diff(pos) < dt, 0]; % Append a 0 at the end to force ending with nothing to merge
iNear = find( diff(x) == 1 ) + 1; % Ex: iNear=3-> merge markers 3 & 4
if x(1)==1, iNear = [1,iNear]; end % Deal with special case of first 

mrk.merge_l = true(1,nMrk);
% Deal with "consecutive too near markers"; select the first one of them. 
for i=1:length(iNear)
    iConsec = iNear(i)+find(x(iNear(i)+1:end)==0,1)-1; % If no consecutive, iConsec=iNear(i)
    iMerge = iNear(i):iConsec+1;
    [mx,imx] = max(feat(iMerge));
    lMerge = false(size(iMerge)); % Logical
    lMerge(imx) = true; % Kill everyone except the one to keep
    mrk.merge_l(iMerge) = lMerge;
end

end

% -------------------------------------------------------------------------
function sds = Representation(sds)
% Represent each class by its centroid, size and expert content

% Inputs
mrk = sds.candidate.merge.marker; % MERGED candidate markers
% % mrk = sds.candidate.marker; % Candidate markers
mrkseg = sds.segment.marker; % All segment markers
classFeat = sds.hierarchy.feature; % Names of features used for class

% Flag expert events
ref = FlagMarkersBuild(sds, sds.flag.expert_l);
mrk = sc_class_marker('Flag', mrk, ref, 'second', 'expert_l');
mrkseg = sc_class_marker('Flag', mrkseg, ref, 'second', 'expert_l');
% Make it logical
mrk.expert_l = mrk.expert_l>0;
mrkseg.expert_l = mrkseg.expert_l>0;
% Number of expert events
Ne = sum(mrk.expert_l);

% Extract matrix of specified features used for classification
featVal = sc_class_marker('Matrix', mrk, classFeat);
Nf = size(featVal,1);

% Compute classes size, centroid and expert-content
Ng = max(mrk.classnum);
centroid = zeros(Nf,Ng);
classSize = zeros(1,Ng);
expert = zeros(1,Ng); % Distribution of expert events across classes
for ii=1:Ng
    idx = mrk.classnum==ii;
    classSize(ii) = sum(idx);
    centroid(:,ii) = mean(featVal(:,idx),2);
    expert(ii) = sum(mrk.expert_l(idx));
end
% Normalize expert content so it sums to 1
if Ne>0, expert = expert/Ne; end

% Renumerate class number in decreasing class size order
[classSize,iSort] = sort(classSize,'descend');
centroid = centroid(:,iSort);
expert = expert(iSort);
numUnordered = mrk.classnum;
for ii=1:Ng
    idx = numUnordered==iSort(ii);
    mrk.classnum(idx) = ii;
end

% Renumerate class number in decreasing expert content order
[expert,iSort] = sort(expert,'descend');
centroid = centroid(:,iSort);
classSize = classSize(iSort);
numUnordered = mrk.classnum;
for ii=1:Ng
    idx = numUnordered==iSort(ii);
    mrk.classnum(idx) = ii;
end

% OUTPUT
% % sds.candidate.marker = mrk;
sds.candidate.merge.marker = mrk; % Merged candidate markers
sds.segment.marker = mrkseg;
sds.class.centroid = centroid;
sds.class.size = classSize;
sds.class.expert = expert;
sds.class.num = mrk.classnum;

end

% -------------------------------------------------------------------------
function sds = Identification(sds)
% Identifies clusters as interpretable classes

% INPUT
mrk = sds.candidate.merge.marker; % MERGED candidate markers
% % mrk = sds.candidate.marker; % Candidate markers
iFeatFreq = strcmpi(sds.hierarchy.feature,'frequency');
iFeatMed = strcmpi(sds.hierarchy.feature,'median');
clustCtr = sds.class.centroid;
expCutoff = sds.class.expertcutoff/100;
clustExpert = sds.class.expert;

% CONSTANTS
% limFreq = median(sds.class.centroid(iFeatFreq,:)); % legacy:13
limFreq = 13;
limMed = 0.02; % A bit in front of Cz line

% Slow (frequency-) and frontal (median+) class number(s)
iClassSlowfront = find(clustCtr(iFeatFreq,:)<limFreq & clustCtr(iFeatMed,:)>=limMed);
% Fast (frequency+) and frontal (median+) class number(s)
iClassFastfront = find(clustCtr(iFeatFreq,:)>=limFreq & clustCtr(iFeatMed,:)>=limMed);
% Slow (frequency-) and posterior (median-) class number(s)
iClassSlowpost  = find(clustCtr(iFeatFreq,:)<limFreq & clustCtr(iFeatMed,:)<limMed);
% Fast (frequency+) and posterior (median-) class number(s)
iClassFastpost  = find(clustCtr(iFeatFreq,:)>=limFreq & clustCtr(iFeatMed,:)<limMed);
% Expert class number(s)
iClassExpert    = 1:find(cumsum(clustExpert)>=expCutoff,1,'first');

% Add class indentification as flags on candidate markers
mrk.slowfront_l = ismember(mrk.classnum,iClassSlowfront);
mrk.fastfront_l = ismember(mrk.classnum,iClassFastfront);
mrk.slowpost_l  = ismember(mrk.classnum,iClassSlowpost);
mrk.fastpost_l  = ismember(mrk.classnum,iClassFastpost);
mrk.expertclass_l  = ismember(mrk.classnum,iClassExpert);

% OUTPUT
% % sds.candidate.marker = mrk;
sds.candidate.merge.marker = mrk; % Merged candidate markers

end

% =========================================================================
% ========== ANALYSIS RESULTS
% =========================================================================

% -------------------------------------------------------------------------
function sFile = ClearEvents(sFile)
% Clear already existing SDS spindles events

fileEvtLabel = {sFile.events.label};
i2rmv = [];
for iEvt=1:numel(fileEvtLabel)
    if ~isempty(strfind(fileEvtLabel{iEvt},'SPIN_CLASS')),     i2rmv = [i2rmv,iEvt]; continue; end
    if ~isempty(strfind(fileEvtLabel{iEvt},'SPIN_EXPERT')),    i2rmv = [i2rmv,iEvt]; continue; end
    if ~isempty(strfind(fileEvtLabel{iEvt},'SPIN_SLOWFRONT')), i2rmv = [i2rmv,iEvt]; continue; end
    if ~isempty(strfind(fileEvtLabel{iEvt},'SPIN_SLOWPOST')),  i2rmv = [i2rmv,iEvt]; continue; end
    if ~isempty(strfind(fileEvtLabel{iEvt},'SPIN_FASTFRONT')), i2rmv = [i2rmv,iEvt]; continue; end
    if ~isempty(strfind(fileEvtLabel{iEvt},'SPIN_FASTPOST')),  i2rmv = [i2rmv,iEvt]; continue; end
end
sFile.events(i2rmv) = [];
end

% -------------------------------------------------------------------------
function sEvents = CreateEvents(sds)
% Return an array of BST event structs to be added to file

% INPUTS
mrk = sds.candidate.merge.marker; % MERGED candidate markers
% % mrk = sds.candidate.marker; % Candidate markers
color  = lines(100);
nClass = max(mrk.classnum);

sEvents = [];
iColor = 0;

% SLOWFRONT SPINDLES
iColor = iColor+1;
mrki = sc_class_marker('Get',mrk,mrk.slowfront_l>0);
% Basic events structure
sEvent = db_template('event');
sEvent.label   = 'SPIN_SLOWFRONT';
sEvent.color   = color(iColor,:);
sEvent.times   = mrki.position;
sEvent.samples = mrki.position_i;
sEvent.epochs  = ones(1, size(sEvent.times,2));
sEvents = [sEvents, sEvent]; % Add to main container

% SLOWPOST SPINDLES
iColor = iColor+1;
mrki = sc_class_marker('Get',mrk,mrk.slowpost_l>0);
% Basic events structure
sEvent = db_template('event');
sEvent.label   = 'SPIN_SLOWPOST';
sEvent.color   = color(iColor,:);
sEvent.times   = mrki.position;
sEvent.samples = mrki.position_i;
sEvent.epochs  = ones(1, size(sEvent.times,2));
sEvents = [sEvents, sEvent]; % Add to main container

% FASTFRONT SPINDLES
iColor = iColor+1;
mrki = sc_class_marker('Get',mrk,mrk.fastfront_l>0);
% Basic events structure
sEvent = db_template('event');
sEvent.label   = 'SPIN_FASTFRONT';
sEvent.color   = color(iColor,:);
sEvent.times   = mrki.position;
sEvent.samples = mrki.position_i;
sEvent.epochs  = ones(1, size(sEvent.times,2));
sEvents = [sEvents, sEvent]; % Add to main container

% FASTPOST SPINDLES
iColor = iColor+1;
mrki = sc_class_marker('Get',mrk,mrk.fastpost_l>0);
% Basic events structure
sEvent = db_template('event');
sEvent.label   = 'SPIN_FASTPOST';
sEvent.color   = color(iColor,:);
sEvent.times   = mrki.position;
sEvent.samples = mrki.position_i;
sEvent.epochs  = ones(1, size(sEvent.times,2));
sEvents = [sEvents, sEvent]; % Add to main container

% EXPERT SLOWFRONT SPINDLES
iColor = iColor+1;
idx = mrk.slowfront_l>0 & mrk.expertclass_l>0;
mrki = sc_class_marker('Get',mrk,idx);
% Basic events structure
sEvent = db_template('event');
sEvent.label   = 'SPIN_SLOWFRONT_EXPERT';
sEvent.color   = color(iColor,:);
sEvent.times   = mrki.position;
sEvent.samples = mrki.position_i;
sEvent.epochs  = ones(1, size(sEvent.times,2));
sEvents = [sEvents, sEvent]; % Add to main container

% EXPERT SLOWPOST SPINDLES
iColor = iColor+1;
idx = mrk.slowpost_l>0 & mrk.expertclass_l>0;
mrki = sc_class_marker('Get',mrk,idx);
% Basic events structure
sEvent = db_template('event');
sEvent.label   = 'SPIN_SLOWPOST_EXPERT';
sEvent.color   = color(iColor,:);
sEvent.times   = mrki.position;
sEvent.samples = mrki.position_i;
sEvent.epochs  = ones(1, size(sEvent.times,2));
sEvents = [sEvents, sEvent]; % Add to main container

% EXPERT FASTFRONT SPINDLES
iColor = iColor+1;
idx = mrk.fastfront_l>0 & mrk.expertclass_l>0;
mrki = sc_class_marker('Get',mrk,idx);
% Basic events structure
sEvent = db_template('event');
sEvent.label   = 'SPIN_FASTFRONT_EXPERT';
sEvent.color   = color(iColor,:);
sEvent.times   = mrki.position;
sEvent.samples = mrki.position_i;
sEvent.epochs  = ones(1, size(sEvent.times,2));
sEvents = [sEvents, sEvent]; % Add to main container

% EXPERT FASTPOST SPINDLES
iColor = iColor+1;
idx = mrk.fastpost_l>0 & mrk.expertclass_l>0;
mrki = sc_class_marker('Get',mrk,idx);
% Basic events structure
sEvent = db_template('event');
sEvent.label   = 'SPIN_FASTPOST_EXPERT';
sEvent.color   = color(iColor,:);
sEvent.times   = mrki.position;
sEvent.samples = mrki.position_i;
sEvent.epochs  = ones(1, size(sEvent.times,2));
sEvents = [sEvents, sEvent]; % Add to main container

% EXPERT SPINDLES
iColor = iColor+1;
idx = mrk.expertclass_l>0;
mrki = sc_class_marker('Get',mrk,idx);
% Basic events structure
sEvent = db_template('event');
sEvent.label   = 'SPIN_EXPERT';
sEvent.color   = color(iColor,:);
sEvent.times   = mrki.position;
sEvent.samples = mrki.position_i;
sEvent.epochs  = ones(1, size(sEvent.times,2));
sEvents = [sEvents, sEvent]; % Add to main container

% Create events for each class
for iClass=1:nClass
    mrki = sc_class_marker('Get',mrk,mrk.classnum==iClass);
    % Basic events structure
    sEvent = db_template('event');
    sEvent.label   = ['SPIN_CLASS',num2str(iClass)];
    sEvent.color   = [0,0,0];
    sEvent.times   = mrki.position;
    sEvent.samples = mrki.position_i;
    sEvent.epochs  = ones(1, size(sEvent.times,2));
    sEvents = [sEvents, sEvent]; % Add to main container
end

end

% =========================================================================
% ========== SHOW RESULTS
% =========================================================================

function ShowMainFigure(sds,figName)

NL = 4; NC = 2;

figure('Name',figName,'position',[839 144 709 803]);
isub=1; % Subplot index
% Show histograms
subplot(NL,NC,isub); isub = isub+1;
ShowHistogram(sds,'segment');
subplot(NL,NC,isub); isub = isub+1;
ShowHistogram(sds,'expert');
subplot(NL,NC,isub); isub = isub+1;
ShowHistogram(sds,'candidate');
subplot(NL,NC,isub); isub = isub+1;
ShowHistogram(sds,'expert-candidate');
% Show resulting centroids in classification features space (optional)
subplot(NL,NC,isub); isub = isub+1;
ShowCentroid(sds);
% Show thresholded dendrogram
subplot(NL,NC,isub); isub = isub+1;
ShowDendrogram(sds);
% Show spatial (2D) histogram
subplot(NL,NC,isub); isub = isub+1;
ShowSpatialHist(sds,'candidate');
% Show individual channels histogram
subplot(NL,NC,isub); isub = isub+1;
ShowChannelHist(sds,'candidate');

end

% -------------------------------------------------------------------------
function ShowCentroid(sds,varargin)
% Any optional parameter is sent to SCATTER.

clustCentroid = sds.class.centroid;
clustSize     = sds.class.size;
featBin       = sds.hierarchy.bin;
featName      = sds.hierarchy.feature;
% % % % clustNum      = sds.hierarchy.class;
% % % % featVal       = sds.hierarchy.vector;
mrk      = sds.candidate.merge.marker; % MERGED Candidate markers
clustNum = mrk.classnum;
featVal  = sc_class_marker('Matrix', mrk, featName);

nClust = max(clustNum);

% Scatter points classes ----------------

clustSizeDisp = 2000*clustSize/max(clustSize); % Scale for display purpose

hold on;

% Show all points
ms = 50;
% if size(featVal,2)>200, ms = 5; else ms = 50; end
scatter(featVal(1,:),featVal(2,:),ms,clustNum,'.');
% % % % % plot(featVal(1,mrk.expert_l>0),featVal(2,mrk.expert_l>0),'*g');

% Colorbar with title
set(get(colorbar,'ylabel'),'string','class #','rotation',-90,'verticalalignment','baseline');

% Show centroids as circle (radius is size)
scatter(clustCentroid(1,:),clustCentroid(2,:),clustSizeDisp,'k',varargin{:});
for ic=1:nClust, 
    txt = text(clustCentroid(1,ic),clustCentroid(2,ic),sprintf('%d (%d)',ic,clustSize(ic)));
    set(txt,'HorizontalAlignment','center');
end
xlim([min(featBin{1}), max(featBin{1})]);
ylim([min(featBin{2}), max(featBin{2})]);
xlabel(featName{1});
ylabel(featName{2});
title('Candidates classification')

hold off

end

% -------------------------------------------------------------------------
function ShowDendrogram(sds)

linktree  = sds.hierarchy.link;
cutoff    = sds.class.cutoff;
criterion = sds.class.criterion;
method    = sds.hierarchy.method;
N = size(linktree,1);

% hD = dendrogram(linktree, N, 'colorthreshold',cutoff); % Colors wont match those in scatter plot
dendrogram(linktree, N);
title('Candidates hierarchical dendrogram')
xlabel('candidates');
ylabel(sprintf('linkage [%s] distance',method));

switch lower(criterion)
    case 'distance'
        hold on;
        plot(0:N,cutoff*ones(1,N+1),'r');
        hold off;
end

end

% -------------------------------------------------------------------------
function ShowHistogram(sds,class)

if nargin<2, class = 'candidate'; end

featBin       = sds.hierarchy.bin;
featName      = sds.hierarchy.feature;
% featHist      = sds.hierarchy.hist;
mrk = sds.candidate.merge.marker; % MERGED candidate markers
% % mrk = sds.candidate.marker; % Candidate markers
mrkall        = sds.segment.marker; % ALL segment markers

switch lower(class)
    case 'segment' % Exclude markers at frequency edges
        mrki = sc_class_marker('Get',mrkall,mrkall.freqedge_l==0);
        axeTitle = 'All events'; 
    case 'candidate'
% %         mrki = sc_class_marker('Get',mrkall,mrkall.candidate_l>0);
        mrki = mrk;
        axeTitle = 'Candidates selected by threshold'; 
    case 'expert'
        mrki = sc_class_marker('Get',mrkall,mrkall.expert_l>0);
        axeTitle = 'Expert events'; 
    case 'expert-candidate'
        mrki = sc_class_marker('Get',mrk,mrk.expertclass_l>0);
        axeTitle = 'Expert supervised class selection'; 
    otherwise
        error('Invalid histogram class selection');
end
% Add number of events in title
N = numel(mrki.position_i);
axeTitle = sprintf('%s (%d)',axeTitle,N);

% Compute 2D histogram of features
featVal = sc_class_marker('Matrix',mrki,featName);
[featHist,C] = hist3(featVal', featBin);

imagesc(featBin{1},featBin{2},featHist');
title(axeTitle)
xlabel(featName{1});
ylabel(featName{2});
axis xy

% Colorbar with title
set(get(colorbar,'ylabel'),'string','histogram count','rotation',-90,'verticalalignment','baseline');

end

% -------------------------------------------------------------------------
function ShowSpatialHist(sds,class)

if nargin<2, class = 'candidate'; end

mrk    = sds.candidate.merge.marker; % MERGED candidate markers
mrkall = sds.segment.marker;         % ALL segment markers

switch lower(class)
    case 'segment' % Exclude markers at frequency edges
        mrki = sc_class_marker('Get',mrkall,mrkall.freqedge_l==0);
        axeTitle = 'All events'; 
    case 'candidate'
% %         mrki = sc_class_marker('Get',mrkall,mrkall.candidate_l>0);
        mrki = mrk;
        axeTitle = 'Candidates selected by threshold'; 
    case 'expert'
        mrki = sc_class_marker('Get',mrkall,mrkall.expert_l>0);
        axeTitle = 'Expert events'; 
    case 'expert-candidate'
        mrki = sc_class_marker('Get',mrk,mrk.expertclass_l>0);
        axeTitle = 'Expert supervised class selection'; 
    otherwise
        error('Invalid histogram class selection');
end
% Add number of events in title
N = numel(mrki.position_i);
axeTitle = sprintf('%s (%d)',axeTitle,N);

% Get selected class spatial features
% Coordinates:
% Row 1: front(+)-back(-)
% Row 2: left(+)-right(-) (weird but yes)
% Row 3: up(+)-down(-)
backfront = mrki.coordinates(1,:);
leftright = -mrki.coordinates(2,:); % Inverse to get left(-)-right(+)

% Compute 2D histogram of features
featVal = [backfront; leftright];
featBin = {-.1:.05:.1, -.1:.05:.1};
[featHist,C] = hist3(featVal', featBin);

imagesc(featBin{2},featBin{1},featHist);
title(axeTitle)
xlabel('left-right position (m)');
ylabel('back-front position (m)');
axis xy

% Colorbar with title
set(get(colorbar,'ylabel'),'string','histogram count','rotation',-90,'verticalalignment','baseline');

end

% -------------------------------------------------------------------------
function ShowChannelHist(sds,class)

if nargin<2, class = 'candidate'; end

mrk    = sds.candidate.merge.marker; % MERGED candidate markers
mrkall = sds.segment.marker;         % ALL segment markers

switch lower(class)
    case 'segment' % Exclude markers at frequency edges
        mrki = sc_class_marker('Get',mrkall,mrkall.freqedge_l==0);
        axeTitle = 'All events'; 
    case 'candidate'
% %         mrki = sc_class_marker('Get',mrkall,mrkall.candidate_l>0);
        mrki = mrk;
        axeTitle = 'Candidates selected by threshold'; 
    case 'expert'
        mrki = sc_class_marker('Get',mrkall,mrkall.expert_l>0);
        axeTitle = 'Expert events'; 
    case 'expert-candidate'
        mrki = sc_class_marker('Get',mrk,mrk.expertclass_l>0);
        axeTitle = 'Expert supervised class selection'; 
    otherwise
        error('Invalid histogram class selection');
end
% Add number of events in title
N = numel(mrki.position_i);
axeTitle = sprintf('%s (%d)',axeTitle,N);

% Get selected class channel info
labels = sds.channel.name;
channel = mrki.channel;
channel_i = cellfun(@(x)find(strcmpi(x,labels)), channel);

% Compute channel histogram
nChan = numel(labels);
[h,b] = hist(channel_i,1:nChan);

% Plot bar histogram
bar(b,h);
set(gca,'xtick',b,'xticklabel',sds.channel.name(b));
rotateXLabels( gca, 90);
title(axeTitle);
ylabel('histogram count');

end

% -------------------------------------------------------------------------
function ExportCSV(sds,filename)

mrk = sds.candidate.merge.marker; % MERGED candidate markers
% % mrk = sds.candidate.marker; % Candidate markers
nMrk = numel(mrk.position);

% Add amplitude gain to make units uV
mrk.amplitude = 10^6*mrk.amplitude;

delimiter = ';';
featName   = {'position','position_i','channel','amplitude','sigmaindex','frequency','median','rightleft','classnum','expert_l','slowfront_l','slowpost_l','fastfront_l','fastpost_l','expertclass_l'};
featFormat = {'%.3f',    '%d',        '%s',     '%.3f',     '%.3f',      '%.3f',     '%.3f',  '%.3f',     '%d',      '%d',      '%d',         '%d',        '%d',         '%d',        '%d'           };

iCell = 3; % Features that are cell arrays of strings 

% Open the file to write in it
fid = fopen(filename,'wb');

Nr = nMrk;
Nc = numel(featName);

% Write the HEADER
iColFormat = ['%s',delimiter];
for iCol=1:Nc
    if iCol==Nc, iColFormat = '%s'; end
    fprintf(fid, iColFormat, featName{iCol});
end
fprintf(fid, '\r\n');

% Write the CSV data
for iRow = 1:Nr % Start at 2 since 1 is header
    for iCol = 1:Nc
        if iCol==Nc, iColFormat = [featFormat{iCol}];
        else         iColFormat = [featFormat{iCol},delimiter];
        end        
        if ismember(iCol,iCell)
            fprintf(fid,iColFormat, mrk.(featName{iCol}){iRow});
        else
            fprintf(fid,iColFormat, mrk.(featName{iCol})(iRow));
        end
    end
    fprintf(fid, '\r\n');
end

fclose(fid);

% % % % % % % Write header
% % % % % % fprintf(fid,   'position; position_i; channel; amplitude; median; frequency; sigmaindex; expert_l; classnum; slowfront_l; slowpost_l; fastfront_l; fastpost_l; expertclass_l');
% % % % % % fprintf(fid, '\r\n');
% % % % % % fprintf
% % % % % % 
% % % % % % % Write the CSV data
% % % % % % for iMrk = 1:nMrk 
% % % % % %     fprintf(fid,'%.3f; %d; ', ...
% % % % % %         );


end

% =========================================================================
% ========== MISC
% =========================================================================

% -------------------------------------------------------------------------
function mrk = FlagMarkersBuild(sds, evtName)
% Create a Marker class struct (format of SC_CLASS_MARKER) from existing
% event in file named EVTNAME. The latter must be extended events,
% otherwise a default 500ms extension is added (to the created markers, not
% the actual events in file). 

FileEvtName = {sds.file.events.label};
mrk.start_i = [];
mrk.end_i = [];
for iEvt=1:numel(evtName)
    idx = find(strcmpi(FileEvtName,evtName{iEvt}));
    if isempty(idx), continue; end
    mrk.start_i = [mrk.start_i, sds.file.events(idx).samples(1,:)];
    if size(sds.file.events(idx).samples,1)==2
        mrk.end_i = [mrk.end_i, sds.file.events(idx).samples(2,:)];
    else % Must add artificial %%%%%%%%%%%%%%%%%
%         mrk.end_i = [mrk.end_i, sds.file.events(idx).samples(1,:)];
        warning('Default duration of 0.5s was used for events [%s]',evtName{iEvt});
        DEFAULT_LENGTH = round(sds.file.prop.sfreq * 0.5);
        mrk.end_i = [mrk.end_i, sds.file.events(idx).samples(1,:)+DEFAULT_LENGTH];
    end
end
mrk.position_i = mrk.start_i;

mrk = sc_class_marker('Remove',mrk, mrk.start_i==0);

end

% -------------------------------------------------------------------------
function EvtNames = ExtractEvtNames(EvtStr)

EvtNames = EvtStr;
if ~isempty(EvtNames)
    EvtNames  = textscan(EvtNames,'%s','delimiter',','); 
    EvtNames = EvtNames{1};
end
    
end

% =========================================================================
% ========== LEGACY
% =========================================================================

% -------------------------------------------------------------------------
function [sds] = FlagMarkers(sds)

mrk = sds.segment.marker;

flagNames = fieldnames(sds.flag);
Nf = numel(flagNames);
for iFlag=1:Nf
    % Create a marker from file event
    ref = FlagMarkersBuild(sds, sds.flag.(flagNames{iFlag}));
    % Flag segments with that marker
    mrk = sc_class_marker('Flag', mrk, ref, 'second', flagNames{iFlag});
end

% Null hypothesis markers flag
if isempty(sds.flag.h0_l)
    mrk.h0_l = true(size(mrk.h0_l));
else
    mrk.h0_l = (mrk.h0_l>0);
end

% Flag to exlude markers outside some markers (at the end of detection)
if isempty(sds.flag.excludeout_l)
    mrk.excludeout_l = false(size(mrk.excludeout_l));
else
    mrk.excludeout_l = ~(mrk.excludeout_l>0);
end

% Flag to reject markers (artefacts) outside some markers
if isempty(sds.flag.rejectout_l)
    mrk.rejectout_l = false(size(mrk.rejectout_l));
else
    mrk.rejectout_l = ~(mrk.rejectout_l>0);
end

% Flag to reject markers (artefacts) inside some markers
mrk.rejectin_l = (mrk.rejectin_l>0);

% Flag events at the edge of frequency band
mrk.freqedge_l = ismember(mrk.frequency, sds.frequency.sigma);

sds.segment.marker = mrk;

% Check if ok
if all(mrk.excludeout_l | mrk.rejectout_l | mrk.rejectin_l | mrk.freqedge_l)
    errMsg = 'All markers are rejected!';
    sds = errMsg;
    return;
end

end

% -------------------------------------------------------------------------
function sds = MergeOld(sds)
% Maybe should be run after all, to merge only markers from same class

if ~ismember(sds.merge.feature, sds.threshold.feature)
    error('Must define a feature that was used for selection'); 
end

candId = sds.segment.marker.candidate_l;
feat  = sds.segment.marker.(sds.merge.feature)(candId);
pos   = sds.segment.marker.position(candId);

candId = candId(candId); % Select only unmerged candidates

% Find those markers that are too near and deal with "consecutive too near
% markers"; select the first one of them. Ex: iMerge=3 -> the 3rd and 4th
% markers need to be merged
x = diff(pos) < sds.merge.interval;
iNear = find( diff(x) == 1 ) + 1;

for i=1:length(iNear)
    iConsec = iNear(i)+find(x(iNear(i)+1:end)==0,1)-1; % If no consecutive, iConsec=iNear(i)
    iMerge = iNear(i):iConsec+1;
    [mx,imx] = max(feat(iMerge));
    lMerge = false(size(iMerge)); % Logical
    lMerge(imx) = true;
    candId(iMerge) = lMerge;
end

sds.segment.marker.candidate_l(sds.segment.marker.candidate_l) = candId;

end

% -------------------------------------------------------------------------
function sds = Topography(sds)

% Get candidate markers
mrk = sds.candidate.marker;
% Channels index to consider for topography
iChan = sds.channel.idx; 
% Duration (samples) of data window for topo computation
duration_i = round(sds.topography.duration*sds.file.prop.sfreq);
% Filtering options
LowPass = sds.frequency.sigma(2);
HighPass = sds.frequency.sigma(1);
FS = sds.file.prop.sfreq;

nMrk = sc_class_marker('Numel',mrk); % Number of candidates
nChan = numel(iChan); % Number of channels
% Samples bounds around markers
SampleBounds = mrk.position_i - round(duration_i/2);
SampleBounds = [SampleBounds; (mrk.position_i + round(duration_i/2))];

% Loop and compute on all markers
bst_progress('start', 'Feature Extraction', 'Computing topographies', 0, nMrk);
mrk.topography = zeros(nChan,nMrk);
mrk.coordinates = zeros(3,nMrk);
for iMrk=1:nMrk
    bst_progress('set', iMrk);
    
    % Get window data
    [F, TimeVector] = in_fread(sds.file, ...
        1, ... % epoch
        SampleBounds(:,iMrk), ...
        iChan);
    % Filter in specific band
    F = detrend(F','linear')';
    F = bst_bandpass(F, FS, HighPass, LowPass);
    % SVD: topographic components (U vectors)
    [U,S,V] = svd(F);
    % Topography is the first component
    mrk.topography(:,iMrk) = U(:,1);
end
% Get location at maximum of topography
% Row 1: front(+)-back(-)
% Row 2: left(+)-right(-) (weird but yes)
% Row 3: up(+)-down(-)
[maxTopo,idxTopo] = max(abs(mrk.topography));
mrk.coordinates = sds.channel.loc(:,idxTopo); 

% Median position (front(+)-back(-))
mrk.median = mrk.coordinates(1,:);

sds.candidate.marker = mrk;

end

% =========================================================================
% ========== JUNK
% =========================================================================

% % % %         % Channel index
% % % %         ChanNames = {sChan.Channel.Name}; % All channel names in file
% % % %         % Get channels names specified
% % % %         SelChanNames = textscan(sProcess.options.sensortypes.Value,'%s','delimiter',',');
% % % %         SelChanNames = SelChanNames{1};
% % % %         % Index of channels specified not present in file
% % % %         iNotFound = find(~ismember(lower(SelChanNames),lower(ChanNames)));
% % % %         if ~isempty(iNotFound)
% % % %             msgstr = 'Following channels are unavailable: ';
% % % %             msgstr = [msgstr,sprintf('\n%s',SelChanNames{iNotFound})];
% % % %             bst_report('Error', sProcess, sInputs(iInput), msgstr);
% % % %             continue;
% % % %         end
% % % %         % Set indexes of selected channels
% % % %         Nc = numel(channelNames);
% % % %         iChannels = zeros(1,Nc);
% % % %         for i=1:Nc
% % % %             iChannels(i) = find(ismember(lower(ChanNames),lower(SelChanNames(i))));
% % % %         end
