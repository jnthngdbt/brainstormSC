
% Ask user to select a Stellate SIG file to test-open
[filename, pathname, filterindex] = uigetfile('*.sig', 'Pick a Stellate Harmonie SIG file');
Filename = [pathname,filename];

iGrpList = 1;
iMtgList = 1;

% open a file and keep it open in MATLAB environment till
% it is closed using mFileClose. Using a interface function
% with a new file nemae will also close the current open file 
% and open the new one in the environment
mFileOpen(Filename);

% close the file in MATLAB environment
mFileClose;

% check the file length
[NumRecs, NumSamps, NumSecs] = mGetFileLength(Filename);

% check the number of recording channels
[NumOfChan] = mGetNumChan(Filename);

% check the base sampling frequency
[TrueSampFreq] = mGetTrueSampFreq(Filename);

% the time of the first recorded sample in the given SIG file.
[RecordingStartTime] = mGetRecStartTime(Filename);

% read the montage list
[MtgList] = mGetMtgList(Filename);

% read the electrode names of the recording montage
[RecChanNames] = mGetRecChanName(Filename);

% extracts the reformatting montage formula
[ChanFormula]= mGetMtgChanFormula(Filename, MtgList{iMtgList});

% read the montage channel labels
[ChanLabel]= mGetChanLabel(Filename, MtgList{iMtgList});

% check the montage channel types
[ChanType]= mGetChanType(Filename, MtgList{iMtgList});

% check the sampling rate of every recording channel 
[SampRate] = mGetChanSampRate(Filename);

% check the event group list
[GrpList] = mGetEvtGrpList(Filename);

% read the number of events in one event group
[NumStatusItems] = mGetNumStatusItemsOfEvt(Filename,GrpList{iGrpList});

% retrieves all the status items and their properties 
[PropHeader, StatusItems] = mGetStatusItems(Filename, GrpList{iGrpList});


% retrieves all the status items start time
[StatusItemTime] = mGetStatusItemTime(Filename,GrpList{iGrpList});

%  retrieves the status items start time of the 2nd event item
[StatusItemTime2] = mGetStatusItemTimeAt(Filename, 3, GrpList{iGrpList});

% the corresponding detecting montage name
[MtgName] = mGetDetMtg(Filename, GrpList{iGrpList});

% read the EEG data and the sampling rate of each channel given offset and
% duration in samples
[EegData, SampRate1] = mGetDataSamp(Filename, 0, 500, MtgList{iMtgList});

% read the EEG data and the sampling rate of each channel given offset and
% duration in seconds
[EegData1, SampRate2] = mGetDataSec(Filename, 0, 2.5, MtgList{iMtgList});

% read the EEG data and the sampling rate of each channel given the start
% time and duration
[EegData2, SampRate3] = mGetDataTime(Filename, '11:02:15.00', 2.5, MtgList{iMtgList});

mFileClose;

disp('The MATLAB SIG File Interface 6.2 was successfully installed in your computer!');
