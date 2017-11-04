%MGETTRUESAMPFREQ.M
%
% This function get base sampleing frequency in the given sig file.
%
% Usage:
%   [TrueSampFreq] = mGetTrueSampFreq(Filename)
%
% Input:
%   Filename-     name of the SIG file. Default extension is '.sig'. 
%                 This file must be in MATLAB search path or current working 
%                 directory.
% Output:
%   TrueSampFreq- the base sampling frequency
%             
%
% See also: mGetChanSampRate
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [TrueSampFreq] = mGetTrueSampFreq(Filename)

NumOfChan = 0;
MFilename = lower(mfilename);
%disp(MFilename);


if MFilename == 'mgettruesampfreq',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 1)
       beep
       disp('mGetTrueSampFreq requires one input arguments.');
       return;
   end
   
   if (nlhs > 1) 
       beep
       disp('mGetTrueSampFreq requires one output argument.');
       return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(Filename)| size(Filename,1)~=1 )
      beep
      disp('mGetTrueSampFreq requires that Sig Filename be one character string.');
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
  
  % call the MEX function
  [TrueSampFreq] =  mSigFileInterFace(MFilename, Filename);
  
% return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
