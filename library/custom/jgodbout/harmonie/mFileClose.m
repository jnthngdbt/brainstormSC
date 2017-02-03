%MFILECLOSE.M
%
% This function is to close a previously opened SIG file.
% 
%
% Usage:
%   mFileClose
%
% 
% See also: mFileOpen
%
% Last modified: Apr. 16, 2003
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mFileClose

NumOfChan = 0;
MFilename = lower(mfilename);
%disp(MFilename);


if MFilename == 'mfileclose',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 0)
      beep
      disp('mFileClose requires no input argument.');
      return;
   end
  
  
 
  % call MEX function
  mSigFileInterFace(MFilename);
  
% return with message indicating invalid function
else
   beep
   disp(['Invalid function name: ' MFileName]);
end
