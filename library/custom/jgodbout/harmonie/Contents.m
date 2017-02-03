% Stellate SIG file interface Toolbox.
% Version 1.0   Stellate Systems Inc.   Jan. 2003
% 
% File Open/Close 
%   mFileOpen                - explicitly open a sig file for other
%                              manipulation
%   mFileClose               - explicitly close a SIG file (opened using
%                              mFileOpen or or other command as follows)
%
% Get general information about the SIG file.
%   mGetFileLength           - get file length in number of records, samples and seconds
%   mGetNumChan              - get the number of channels of the recording montage
%   mGetTrueSampFreq         - get the base sampling frequency
%   mGetRecStartTime         - get the start time of the recording
%   mGetMtgList              - get names of all the montages associated with the file 
%   mGetRecChanName          - get recording montage electrode names
%   mGetMtgChanFormula       - get the reformating formula of a specific reformating montage 
%   mGetChanLabel            - get the channel labels of a specific montage
%   mGetChanType             - get the trace type of each channel in a specific montage
%   mGetChanSampRate         - get the sampling rate of each channel
% 
% Get event related information.
%   mGetEvtGrpList           - retrieve event group list
%   mGetNumStatusItemsOfEvt  - retrieve the number of events in a specific event group
%   mGetStatusItems          - retrieve all the status items and their properties of a specific event
%   mGetStatusItemTime       - retrieve start time of all the status items of a specific event
%   mGetStatusItemTimeAt     - retrieve start time of the status item at one specific index
%   mGetDetMtg               - retrieve the montage name used to detect a specific event group
%   mGetDetChan              - retrieve names of detection channels
%
% Read data from SIG file.
%   mGetDataSamp             - Read data block given the offset and the duration in number of samples
%   mGetDataSec              - Read data block given the offset and the duration in number of seconds
%   mGetDataTime             - Read data block given the duration in number of seconds and the start time 
%
% Change/add event items.
%   mChangeStatusItemEvtName - change the name of specific events in one event group
%   mMarkEvents              - mark new event item to the SIG file
%   mModifyEventItem         - change event property value
%
% Create new SIG file.
%   mCreateModifiedSigFile   - create new SIG file based on an existing SIG file
%

                                                 
                      
                        


