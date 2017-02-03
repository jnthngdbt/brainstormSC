%MGETNUMCHAN.M
%
% This function get the number of channels of 
% data recorded in the given SIG file.
%
% Usage:
%   [NumOfChan] = mGetNumChan(Filename)
%
% Input:
%   Filename-  name of the SIG file. Default extension is '.sig'. 
%              This file must be in MATLAB search path or current working 
%              directory.
% Output:
%   NumOfChan- number of channels of the recording montage
%             
%
% See also: mGetRecChanName, mGetTrueSampFreq, mGetMtgList, mGetMtgChanFormula
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NumOfChan] = mGetNumChan(Filename)

NumOfChan = 0;
MFilename = lower(mfilename);
%disp(MFilename);


if MFilename == 'mgetnumchan',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 1)
      beep
      disp('mGetNumChan requires one input arguments.');
      return;
   end
	if (nlhs > 1) 
      beep
      disp('mGetNumChan requires one output argument.');
      return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(Filename) | size(Filename,1)~=1 )
     beep
     disp('mGetNumChan requires that Sig Filename be one character string.');
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
  
  % call MEX function
  [NumOfChan] =  mSigFileInterFace(MFilename, Filename);
  
% return with message indicating invalid function
else
   beep
   disp(['Invalid function name: ' MFileName]);
end
