function varargout = process_sc_sleep_cycle( varargin )
% process_sc_sleep_cycle: Create sleep cycle markers from sleep scoring
% markers. 
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'SLEEP > Sleep cycles from REM periods';
    sProcess.FileTag     = '';
    sProcess.Description = '<HTML>Sleep cycles';
    sProcess.Category    = 'custom';
    sProcess.Index       = sc_bst_process_index(2100);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % 
    sProcess.options.REM.Comment = 'REM events name: ';
    sProcess.options.REM.Type    = 'text';
    sProcess.options.REM.Value   = 'REM';
    %
    sProcess.options.MRS.Comment = 'Minimum time of REM state to start a REM period (MRS): ';
    sProcess.options.MRS.Type    = 'value';
    sProcess.options.MRS.Value   = {0,'min',3};
    %
    sProcess.options.MRE.Comment = 'Minimum time of not-REM state to end a REM period (MRE): ';
    sProcess.options.MRE.Type    = 'value';
    sProcess.options.MRE.Value   = {15,'min',3};
    % === 
    sProcess.options.label2.Comment = [ ...
        '<HTML><BR>' ...
'<BR> A REM period starts when the MRS minutes following a REM marker (NOT in a' ...
'<BR> current REM period) is also marked as REM. A REM period ends when the MRE' ...
'<BR> minutes following a REM marker (in a current REM period) do NOT contain' ...
'<BR> REM markers.' ...
        ];
    sProcess.options.label2.Type    = 'label';
% '<BR> A REM period starts MRS minutes after the first REM marker appearance' ...
% '<BR> since last REM period. A REM period ends MRE minutes after the end of ' ...
% '<BR> the last REM marker seen. A new cycle begins when a REM period ends.' ...

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    
    bst_progress('start', 'Creating sleep cycles', 'Creating sleep cycles');
    
    % Get parameters
    MRS = sProcess.options.MRS.Value{1};
    MRE = sProcess.options.MRE.Value{1};
    REMName = sProcess.options.REM.Value;
    DT = 10; % Time resolution (s)
    
    for iFile = 1:length(sInputs)

        % ===== GET DATA =====
        % Load the raw file descriptor
        DataMat = in_bst_data(sInputs(iFile).FileName);
        sFile = DataMat.F;
        
        % ===== COMPUTE SLEEP CYCLES FROM REM MARKERS
        % Get REM events
        iREMEvt = find(strcmpi(REMName,{sFile.events.label}));
        if isempty(iREMEvt)
            strMsg = sprintf('Could not find REM events. No events named [%s]',REMName);
            bst_report('Error', sProcess, sInputs(iFile), strMsg);
            continue
        end
        REMBounds = sFile.events(iREMEvt).times;
        % Get cycles time bounds
        CyclesBounds = SleepCyclesFromREM(REMBounds, MRS, MRE, DT, sFile.prop.times(2));
        % Create cycles events
        sEvents = db_template('event');
        sEvents.label   = 'SleepCycles';
        sEvents.color   = [0,0,0];
        sEvents.times   = CyclesBounds;
        sEvents.samples = round(CyclesBounds*sFile.prop.sfreq)+1;
        sEvents.epochs  = ones(1, size(sEvents.times,2));
        % Create 1 event for each cycle
        for iCycle = 1:size(CyclesBounds,2)
            sEvent = db_template('event');
            sEvent.label   = ['SleepCycle',num2str(iCycle)];
            sEvent.color   = [0,0,0];
            sEvent.times   = CyclesBounds(:,iCycle);
            sEvent.samples = round(sEvent.times*sFile.prop.sfreq)+1;
            sEvent.epochs  = 1;
            sEvents = [sEvents, sEvent];
        end

        % ===== SAVE RESULT =====
        % Delete previous cycle events
        fileEvtLabel = {sFile.events.label};
        i2rmv = [];
        for iEvt=1:numel(fileEvtLabel)
            if numel(fileEvtLabel{iEvt})<10, continue; end
            if strcmpi(fileEvtLabel{iEvt}(1:10),'SleepCycle'), i2rmv = [i2rmv,iEvt]; continue; end
        end
        sFile.events(i2rmv) = [];
        % Only save changes if something was detected
        if ~isempty(sEvents)
            % Update events structure
            if ~isfield(sFile, 'events') || isempty(sFile.events)
                sFile.events = sEvents;
            else
                sFile.events = [sFile.events, sEvents];
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
    OutputFiles = {sInputs.FileName};
            
end

% =========================================================================
function CyclesBounds = SleepCyclesFromREM(REMBounds, MRS, MRE, DT, MaxBound)
% A REM period starts when the MRS minutes following a REM marker (NOT in a
% current REM period) is also marked as REM. A REM period ends when the MRE
% minutes following a REM marker (in a current REM period) do NOT contain
% REM markers.
% 
% REM markers are represented by REMBOUNDS, a 2xN matrix (for N REM
% markers) where the first line is the start times (s) of all REM markers
% and the second line is the end times (s) of all REM markers.
% 
% CYCLESBOUNDS is a 2xM matrix, where M is the number of sleep cycles.
% First and second lines are start and end times (s) respectively.
% 
% (optional) DT is the desired time sampling interval (time resolution)
% 
% Default values: MRS=0min, MRE=15min, DT=10s
% 
% Jonathan Godbout, 2013

% Default inputs
if nargin<2, MRS = 0; end  % minutes
if nargin<3, MRE = 15; end % minutes
if nargin<4, DT = 10; end  % seconds
if nargin<5, MaxBound = max(REMBounds(:)); end
FS = 1/DT;

% Time vector, long enough to include all REM period markers
t = 0:DT:MaxBound;
% Create REM state vector (REM-only hypnogram)
REMHypno = false(1,numel(t));
for iREM = 1:size(REMBounds,2)
    REMHypno = REMHypno | (t>=REMBounds(1,iREM) & t<=REMBounds(2,iREM));
end

% Find cycle ends
REMPeriods = FlipFlop(REMHypno, 0, MRS*60*FS, 0.5, 0, MRE*60*FS, 0.5);
iCycleEnd = t(diff(REMPeriods)==-1); % (s)

% Create cycles bounds
CyclesBounds = [0; iCycleEnd(1)]; % (s)
for i=1:numel(iCycleEnd)-1
    NewBounds = [iCycleEnd(i); iCycleEnd(i+1)];
    CyclesBounds = [CyclesBounds, NewBounds];
end
NewBounds = [iCycleEnd(end); MaxBound];
CyclesBounds = [CyclesBounds, NewBounds];

end

% =========================================================================
function y = FlipFlop(x,ub,ua,ut,db,da,dt)
% X is a state signal of [1xN] samples. UB and UA is the number of samples
% before and after sample n (respectively) where X[n]>=UT for Y[n] to be
% set to 1 and keep this state until Y is set to 0 (in other words, the
% condition for an "up" state). Similarly, DB and DA is the number of samp.
% before and after sample n (respectively) where X[n]<=DT for Y[n] to be
% set to 0 and keep this state until Y is set to 1 (in other words, the
% condition for a "down" state). UT and DT are up and down thresholds on X
% respectively.

[Nc, Ns] = size(x);
if Ns==1
    x = x';
    [Nc, Ns] = size(x);
end

y = zeros(size(x));
for ii=1:Nc
    y(ii,:) = flipflop_1d(x(ii,:),ub,ua,ut,db,da,dt);
end


end

function y = flipflop_1d(x,ub,ua,ut,db,da,dt)

Ns = length(x);

edge_left = max([ub db]);
edge_right = max([ua da]);

y = zeros(size(x));
for ii=edge_left+1:Ns-edge_right
    if     all( x(ii-ub:ii+ua) >= ut )
        y(ii) = 1;
    elseif all( x(ii-db:ii+da) <= dt )
        y(ii) = 0;
    else
        if ii>1
            y(ii) = y(ii-1);
        end
    end
end
y(ii:end) = y(ii); % Extend last state to the end
y(1:edge_left+1) = y(edge_left+1); % Extend first state to the start

end








% % % % % % if ~slp.isfeature('sleep_i')
% % % % % %     methods.c_events.sleep.add_feature_sleep_i(slp);
% % % % % % end
% % % % % % 
% % % % % % slpHypno = slp.feature.sleep_i;
% % % % % % slpHypno(slpHypno>5) = 0;
% % % % % % 
% % % % % % % Sleep scoring interval (scoring window size) and equivalent sampling freq
% % % % % % slpDT = mode(diff(slp.feature.position));
% % % % % % slpFS = 1/slpDT;
% % % % % % 
% % % % % % % Find cycle ends
% % % % % % t = signal.filter.flipflop(slpHypno, 0, MRS*60*slpFS, 4.5, 0, MRE*60*slpFS, 4);
% % % % % % iCycleEnd = find(diff(t)==-1);
% % % % % % 
% % % % % % % Add CYCLE_I feature to sleep events
% % % % % % slpCycle_i = ones(1,slp.N);
% % % % % % for i=1:numel(iCycleEnd)-1
% % % % % %     slpCycle_i(iCycleEnd(i)+1:iCycleEnd(i+1)) = i+1;
% % % % % % end
% % % % % % slpCycle_i(iCycleEnd(end)+1:end) = i+1+1;
% % % % % % 
% % % % % % slp.add_feature('cycle_i',slpCycle_i);
% % % % % % 
% % % % % % slp.add_history(input);
% % % % % % 
% % % % % % % e.apply_feature(slp,'cycle_i');