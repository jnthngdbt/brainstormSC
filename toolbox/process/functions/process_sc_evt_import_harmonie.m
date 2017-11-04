function varargout = process_sc_evt_import_harmonie( varargin )
% process_sc_evt_import_harmonie: Import Harmonie events. To do so, convert
% *.sts file into Brainproducts *.vmrk files and uses Brainstorm buit-in
% function to import this supported file.
%
% External call:
% 
%   process_sc_evt_import_harmonie('ConvertHarmonie2Brainproducts');
%   process_sc_evt_import_harmonie('ConvertHarmonie2Brainproducts',FILES);
%   process_sc_evt_import_harmonie('ConvertHarmonie2Brainproducts',FILES,NEWFS);
%   process_sc_evt_import_harmonie('ConvertHarmonie2Brainproducts',..., 'addchannelinfo');
% 
% where FILES is a CELL array of strings (Stellate file names) or a string
% (single Stellate file). With no input, will ask user to select the files
% to convert. The 'addchannelinfo' flag will add the channel name of the
% markers in their "description" field (only if the marker is on a specific
% channel; "description" wont be changed if the marker is on all channels).
% 
% NEWFS is a vector of new sampling frequencies. It can be unique (same new
% frequency will be applied to all input FILES), otherwise it must match
% the number of FILES. Example, if the markers to import have been marked 
% in a downsampled version of the original file (whose original
% sampling frequency is NEWFS), it will create an additional *.NEWFS.vmrk
% marker file whos markers' time position (in samples) will correspond to a
% sampling frequency NEWFS. If NEWFS is the same as the Stellate files to
% convert, this additional *.NEWFS.vmrk wont be created. 
% 
% Author: Jonathan Godbout, 2013 

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'IMPORT > EVENTS > Stellate Harmonie (only under Windows and Matlab 32bits)';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(65);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'raw'};
    sProcess.OutputTypes = {'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 0;
    sProcess.isSeparator = 0;
    % Option: File to import
    sProcess.options.datafile.Comment = 'You will be asked to select Stellate files';
    sProcess.options.datafile.Type    = 'label';
    sProcess.options.datafile.Value   = [];
    % Separator
    sProcess.options.separator1.Type = 'separator';
    sProcess.options.separator1.Comment = ' ';
    %
    sProcess.options.label0.Comment = [ '<HTML>',...
    'Processing comments:'...
    '<BR> - It will first convert *.STS file into *.VMRK files, then find'...
    '<BR>   the right corresponding raw file by name matching.', ...
    '<BR> - If sampling frequencies of matching Harmonie file and input raw'...
    '<BR>   file are not the same, another *.VMRK file is created with resampled', ...
    '<BR>   markers to correctly register with input raw file.', ...
    ];
    sProcess.options.label0.Type    = 'label';
    % Separator
    sProcess.options.separator2.Type = 'separator';
    sProcess.options.separator2.Comment = ' ';
    % Channel
    sProcess.options.addchaninfo.Comment = 'Add channel info in event names';
    sProcess.options.addchaninfo.Type    = 'checkbox';
    sProcess.options.addchaninfo.Value   = 0;
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>

% Get raw files names
RawFileNames = cell(1,numel(sInputs));
for i=1:numel(sInputs)
    FileMat = in_bst_data(sInputs(i).FileName);
    % Raw filename parts
    [rawpath,rawname,rawext] = fileparts(FileMat.F.filename);
    RawFileNames{i} = rawname;
end

% Get Harmonie file names
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
vmrkFiles = cell(Nf,1);
bst_progress('start', 'Importing Stellate Events', 'Processing...', 0, Nf);
for iInput=1:Nf
    bst_progress('set', iInput);
    
    % File name infos
    [p,subjectname,ext] = fileparts(datafiles{iInput});
    
    % Find match in list of input raw files with current events file 
    iMatch = find(cellfun(@(x)~isempty(strfind(upper(datafiles{iInput}),upper(x))),RawFileNames));
    if isempty(iMatch)
        Fs = 0; % No match found; if Fs==0, convert markers w/o resampling
    else
        % Matching raw data descriptor and sampling frequency
        DataMat = in_bst_data(sInputs(iMatch(1)).FileName);
        Fs = DataMat.F.prop.sfreq;
    end
                
    try
        if sProcess.options.addchaninfo.Value>0 % Add channel in event name
            vmrkFiles(iInput) = ConvertHarmonie2Brainproducts(datafiles{iInput},Fs,'addchannelinfo');      
        else
            vmrkFiles(iInput) = ConvertHarmonie2Brainproducts(datafiles{iInput},Fs);      
        end
    catch % ME not supported in 2006b
        strMsg = sprintf('There was a problem converting file [%s]: \n %s',[subjectname,ext]); % ,getReport(ME)
        bst_report('Warning', sProcess, sInputs, strMsg)
    end
end

% Process: IMPORT > EVENTS > From file (name-matched raw-event files)
sFiles = bst_process(...
    'CallProcess', 'process_sc_evt_import_auto', ...
    {sInputs.FileName}, [], ...
    'evtfile', {vmrkFiles, 'BRAINAMP'});

OutputFiles = {sInputs.FileName};

end

% ========================================================================
function vmrkFiles = ConvertHarmonie2Brainproducts(varargin)
% Convert Stellate configuration file (*.sts) to Brainproduct format 
% (*.vhdr,*.vmrk).
% 
%   ConvertHarmonie2Brainproducts;
%   ConvertHarmonie2Brainproducts(FILES);
%   ConvertHarmonie2Brainproducts(FILES,NEWFS);
%   ConvertHarmonie2Brainproducts(..., 'addchannelinfo');
% 
% where FILES is a CELL array of strings (Stellate file names) or a string
% (single Stellate file). With no input, will ask user to select the files
% to convert. The 'addchannelinfo' flag will add the channel name of the
% markers in their "description" field (only if the marker is on a specific
% channel; "description" wont be changed if the marker is on all channels).
% 
% NEWFS is a vector of new sampling frequencies. It can be unique (same new
% frequency will be applied to all input FILES), otherwise it must match
% the number of FILES. Example, if the markers to import have been marked 
% in a downsampled version of the original file (whose original
% sampling frequency is NEWFS), it will create an additional *.NEWFS.vmrk
% marker file whos markers' time position (in samples) will correspond to a
% sampling frequency NEWFS. If NEWFS is the same as the Stellate files to
% convert, this additional *.NEWFS.vmrk wont be created. 

% Option flag: add channel info in name
iAddChannel = find(strcmpi(varargin,'addchannelinfo'));
if ~isempty(iAddChannel)
    varargin(iAddChannel) = [];
    iAddChannel = true;
else
    iAddChannel = false;
end

% Optional new sample frequency
iNewFs = find(cellfun(@isnumeric,varargin));
if ~isempty(iNewFs)
    NewFs = varargin{iNewFs};
    varargin(iNewFs) = [];
    isResampled = true;
    % Cancel flag if empty or null
    if isempty(NewFs)
        isResampled = false;
    elseif numel(NewFs)==1 && NewFs(1)==0
        isResampled = false;
    end
else
    NewFs = [];
    isResampled = false;
end

[fname,pname] = process_sc_import_harmonie('get_files',varargin{:});
if isempty(fname), return; end

Nf = numel(fname);

% Check if number of new Fs match number of files
if isResampled
    if numel(NewFs)==1
        NewFs = NewFs*ones(1,Nf);
    elseif Nf~=numel(NewFs)
        error('Number of new sampling frequencies must be unique (applied to all files) or match number of input files');
    end
else
    NewFs = zeros(1,Nf);
end

% Read Stellate Harmonie and create Brainproducts files
vmrkFiles = cell(Nf,1);
for ii=1:Nf
    fprintf('%.5f %%\n',100*ii/Nf);
    [hdr,mrk] = process_sc_import_harmonie('OpenHarmonie',[pname,fname{ii}]);
    if isempty(hdr), continue; end % An error was catched
    % Write Brainproducts *.vhdr file
    hdr = process_sc_import_harmonie('write_vhdr',pname,fname{ii},hdr);
    % Write Brainproducts *.vmrk file
    process_sc_import_harmonie('write_vmrk',pname,fname{ii},hdr,mrk,iAddChannel);
    
    % Only create additional resampled marker file if sampling changed
    if isResampled && (round(hdr.Fs)~=round(NewFs(ii)))
        mrk = ResampleMarkers(mrk, hdr.Fs, NewFs(ii));
        % Write Brainproducts *.vmrk file of resampled markers
        fnameResample = sprintf('%s.%dHz',fname{ii},round(NewFs(ii)));
        process_sc_import_harmonie('write_vmrk',pname,fnameResample,hdr,mrk,iAddChannel);
        vmrkFiles{ii} = [pname,fnameResample,'.vmrk'];
    else
        vmrkFiles{ii} = [pname,fname{ii},'.vmrk'];
    end
end

end

% -------------------------------------------------------------------------
function mrk = ResampleMarkers(mrk,fsIn,fsOut)

ResampleFactor = fsOut/fsIn;
for i=1:numel(mrk)
    % Deal with positions at first sample
    mrk(i).feature.position_i(mrk(i).feature.position_i==0) = 1;
    iPos1 = mrk(i).feature.position_i==1;
    % Deal with durations=0
    mrk(i).feature.duration_i(mrk(i).feature.duration_i==1) = 0;
    iDur0 = mrk(i).feature.duration_i==0;
    
    % Resample
    mrk(i).feature.position_i = round(mrk(i).feature.position_i * ResampleFactor);
    mrk(i).feature.duration_i = round(mrk(i).feature.duration_i * ResampleFactor);
    
    % Deal again with positions at first sample
    mrk(i).feature.position_i(iPos1) = 1; % First sample is 1
    % Deal again with durations=0
    mrk(i).feature.duration_i(iDur0) = 1; % Null duration is 1 (not 0)
end

end


% =========================================================================
% LEGACY
% Before dealing with changed sampled frequency
% Also, before overwriting files by default

% % for iInput=1:Nf
% %     bst_progress('set', iInput);
% %     % File name infos
% %     [p,subjectname,ext] = fileparts(datafiles{iInput});
% %     
% %     vmrkFiles{iInput} = [datafiles{iInput},'.vmrk'];
% %     
% %     try
% %         % Skip conversion if already done (converted file exists)
% %         if exist(vmrkFiles{iInput},'file')>0
% %             strMsg = sprintf('Stellate file [%s] is already converted to Brainproduct.',[subjectname,ext]);
% %             bst_report('Warning', sProcess, sInputs, strMsg)
% %         else
% %             if sProcess.options.addchaninfo.Value>0 % Add channel in event name
% %                 vmrkFiles(iInput) = ConvertHarmonie2Brainproducts(datafiles{iInput},'addchannelinfo');      
% %             else
% %                 vmrkFiles(iInput) = ConvertHarmonie2Brainproducts(datafiles{iInput});      
% %             end
% %         end
% %     catch ME
% %         strMsg = sprintf('There was a problem converting file [%s]: \n %s',[subjectname,ext],getReport(ME));
% %         bst_report('Warning', sProcess, sInputs, strMsg)
% %     end
% % end
% % 
% % % Process: IMPORT > EVENTS > From file (name-matched raw-event files)
% % sFiles = bst_process(...
% %     'CallProcess', 'process_sc_evt_import_auto', ...
% %     {sInputs.FileName}, [], ...
% %     'evtfile', {vmrkFiles, 'BRAINAMP'});
% % 
% % OutputFiles = {sInputs.FileName};
% % 
% % end
% =========================================================================
