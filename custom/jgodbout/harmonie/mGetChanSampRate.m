%MGETCHANSAMPRATE.M
%
% Given the file name, this function returns the sampling rate of
% every channel in the recording montage.
%
% Usage:
%   [SampRate] = mGetChanSampRate(Filename)
% 
% Input:
%   Filename-  name of the SIG file. Default extension is '.sig'. 
%              This file must be in MATLAB search path or current working 
%              directory.
% Output:
%   SampRate-  sampling rate of each channel [numChannels x 1]
%             
%
% See also: mGetTrueSampFreq
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [SampRate] = mGetChanSampRate(Filename)

SampRate=[];
MFilename = lower(mfilename);
%disp(MFilename);


if MFilename == 'mgetchansamprate',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 1)
      beep
      disp('mGetChanSampRate requires one input arguments.');
      return;
   end
	if (nlhs > 1) 
      beep
      disp('mGetChanSampRate requires 1 output argument.');
      return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(Filename) | size(Filename,1)~=1 )
     beep
     disp('mGetChanSampRate requires that Sig Filename be one character string.');
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
  
  % call MEX function
  [SampRate] =  mSigFileInterFace(MFilename, Filename);
  
  SampRate = SampRate';
  
% return with message indicating invalid function
else
   beep
   disp(['Invalid function name: ' MFileName]);
end
