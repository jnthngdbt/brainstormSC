%MGETDATATIME.M
%
% This function get the EEG data given the start time of 
% the data block and the duration of the data (in seconds)
% in the given SIG file. The sampling rate of each channel
% is also returned as the second output. Default montage
% is the recording montage.
%
% Note: the polarity of the EEG data is the same as the actual recording value.
%
% Usage:
%   [EegData, SampRate] = mGetDataTime(Filename, Time, DurationSec, MtgName)
% 
% Input:
%   Filename-   name of the SIG file. Default extension is '.sig'. 
%               This file must be in MATLAB search path or current working 
%               directory.
%   Time-       start time string in the form of 'hr:mm:ss.milliseconds'
%   DurationSec-duration in seconds
%   Mtgname-    Montage name. Default is recording montage. It is case sensitive.
%
% Output:
%   EegData-    EEG data required [numSamples x numChannels]. If sampling rate
%               is differnt from channel to channel, zeros will be padded to 
%               the channel of lower sampling rate.
%   SampRate-- sampling rate of each channel [1 x numChannels]
%
% See also: mGetRecStartTime, mGetMtgList, mGetDataSamp, mGetDataSec
%
% Last modified: Jan. 24, 2003
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EegData, SampRate] = mGetDataTime(Filename, Time, DurationSec, MtgName)

EegData = [];
Samprate=[];
MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mgetdatatime',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs <3 | nrhs >4)
       beep
       disp('mGetDataTime requires three or four input arguments.');
       disp('USAGE: [EegData] = mGetDataSec(Filename, OffsetSec, DurationSec, MtgName)');
       return;
   end
   
   if (nlhs > 2) 
       beep
       disp('mGetDataTime requires two output arguments.');
       disp('USAGE: [EegData] = mGetDataSec(Filename, OffsetSec, DurationSec, MtgName)');
       return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(Filename)| size(Filename,1)~=1 )
      beep
      disp('mGetDataTime requires that Sig Filename be a character string.');
      return;
  end
  
  % Check the extension (default extension is '.sig')
  [P N E]=fileparts(Filename);
  if isempty(E)
      Filename=[Filename '.sig'];
  end
  
  % chack the existence of the file
  if exist(Filename,'file')==0
     beep
     disp(['File ' Filename ' does NOT exist!']); 
     return;
  end
    
  % check other data type
  if (~ischar(Time))
      beep
      disp('mGetDataTime requires the Time argument to be a string');
      return;
  end
  if (~isnumeric(DurationSec))
      beep
      disp('mGetDataTime requires the DurationSec argument to be numeric.');
      return;
  end
  
  if (DurationSec<=0)
      return;
  end
  
  % check the montage
  numMtg = -1;
  if nrhs ==4
    if (~ischar(MtgName))
        beep
        disp('mGetDataTime requires that MtgName be one character string.');
        return;
    end
    
    mtgList = mGetMtgList(Filename);
    for i=1:size(mtgList,1)
        if strcmp(MtgName, mtgList(i,1))
            numMtg = i-1;
        end
    end
    
    if numMtg<0
        beep
        disp(['The montage ' MtgName ' does NOT exist.']);
        return;
    end 
else
    numMtg = 0;
end
  
  
  % Call the MEX function
  [EegData, SampRate] =  mSigFileInterFace(MFilename, Filename, Time, DurationSec, numMtg);
     
  % test to see whether data was retrieved
  if (~norm(EegData))
      EegData = [];
  end
   
  % return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
