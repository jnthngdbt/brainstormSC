%MFILEOPEN.M
%
% This function is to open a SIG file for subsequent reading.
% 
%
% Usage:
%   flag = mFileOpen(Filename)
%
% Input:
%   Filename-  name of the SIG file. Default extension is '.sig'. 
%              This file must be in MATLAB search path or current working 
%              directory.
% Output:
%   flag-      flag indicating if the file is opened successfully 
%              (1: succeed, 0: fail)
%                           
% See also: mFileClose
%
% Last modified: Jan. 10, 2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function flag = mFileOpen(Filename)

NumOfChan = 0;
MFilename = lower(mfilename);
%disp(MFilename);


if MFilename == 'mfileopen',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 1)
      beep
      disp('mFileOpen requires one input argument.');
      return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(Filename) | size(Filename,1)~=1 )
     beep
     disp('mFileOpen requires that Sig Filename be one character string.');
     return;
  end
  
  % Check the extension (default extension is '.sig')
  [P N E]=fileparts(Filename);
  if isempty(E)
      Filename=[Filename '.sig'];
  end
  
  % check the existence of the file
  if exist(Filename,'file')==0
     beep
     disp(['File ' Filename ' does NOT exist!']); 
     return;
  end
%   Filename = which(Filename);
  
  % call MEX function
  flag = mSigFileInterFace(MFilename, Filename);
  
% return with message indicating invalid function
else
   beep
   disp(['Invalid function name: ' MFileName]);
end
