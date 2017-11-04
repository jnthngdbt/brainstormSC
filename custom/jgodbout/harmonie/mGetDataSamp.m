%MGETDATASAMP.M
%
% Given the file name, montage name, the offset from the start and 
% the duration (both in number of samples), this function returns
% the EEG data and the sampling rate of each channel. Default montage
% is the recording montage. 
%
% Note: the polarity of the EEG data is the same as the actual recording value.
%
% Usage:
%   [EegData, SampRate] = mGetDataSamp(Filename, OffsetSamp, DurationSamp, MtgName)
%
% Input:
%   Filename-     name of the SIG file. Default extension is '.sig'. 
%                 This file must be in MATLAB search path or current working
%                 directory.
%   OffsetSamp-   offset in sample
%   DurationSamp- duration in sample
%   Mtgname -     Montage name. Default is recording montage. It is case sensitive.
% Output:
%   EegData -     EEG data required [numSamples x numChannels]. If sampling 
%                 rate is differnt from channel to channel, zeros will be 
%                 padded to the channel of lower sampling rate.
%   SampRate-     sampling rate of each channel [1 x numChannels]
%
% See also: mGetMtgList, mGetDataSec, mGetDataTime.
%
% Last modified: Jan 23, 2003
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EegData, SampRate] = mGetDataSamp(Filename, OffsetSamp, DurationSamp, MtgName)


EegData = [];
SampRate = [];
MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mgetdatasamp',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs > 4 |nrhs <3 )
       beep
       disp('mGetDataSamp requires three or four input arguments.');
       disp('USAGE: [EegData, SampRate] = mGetDataSamp(Filename, OffsetSamp, DurationSamp, MtgName)');
       return;
   end
   if (nlhs > 2)
       beep
       disp('mGetDataSamp requires two output arguments.');
       disp('USAGE: [EegData, SampRate] = mGetDataSamp(Filename, StartSamp, DurationSamp, MtgName)');
       return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(Filename)| size(Filename,1)~=1)
      beep
      disp('mGetDataSamp requires that Sig Filename be one character string.');
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
%   Filename = which(Filename);
  
  % check the argument data type
  if (~isnumeric(OffsetSamp))
      beep
      disp('mGetDataSamp requires the Offset argument to be numeric in samples');
      return;
  end
  
  if (~isnumeric(DurationSamp))
      beep
      disp('mGetDataSamp requires the Duration argument to be numeric in samples');
      return;
  end
  
  if (OffsetSamp<0)
      DurationSamp=DurationSamp+OffsetSamp;
      OffsetSamp=0;
  end
  if (DurationSamp<=0)
     return;
  end 
  
  % check the montage
  numMtg = -1;
  if nrhs ==4
    if (~ischar(MtgName))
        beep
        disp('mGetDataSamp requires that MtgName be one character string.');
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
  [EegData, SampRate] = mSigFileInterFace(MFilename, Filename, OffsetSamp, DurationSamp, numMtg);
  
  % test to see whether data was retrieved
  if (~norm(EegData))
     EegData = [];
  end
   
% return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
