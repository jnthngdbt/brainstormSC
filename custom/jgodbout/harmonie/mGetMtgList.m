%MGETMTGLIST.M
%
% Given the file name, this function returns the montage list.
%
% Usage:
%   [MtgList] = mGetMtgList(Filename)
%
% Input:
%   Filename-  name of the SIG file. Default extension is '.sig'. 
%              This file must be in MATLAB search path or current working 
%              directory.
% Output:
%   MtgList -  [numMtgs x 1] cell array of montage names. The first one is 
%              the recording montage while others are reformating montages.
%             
%
% See also: mGetDetMtg, mGetMtgChanFormula, mGetChanLabel
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MtgList] = mGetMtgList(Filename)

NumOfChan = 0;
MFilename = lower(mfilename);
%disp(MFilename);


if MFilename == 'mgetmtglist',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 1)
      beep
      disp('mGetMtgList requires one input arguments.');
      return;
   end
	if (nlhs > 1) 
      beep
      disp('mGetMtgList requires 1 output argument.');
      return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(Filename) | size(Filename,1)~=1 )
     beep
     disp('mGetMtgList requires that Sig Filename be one character string.');
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
  [MtgList] =  mSigFileInterFace(MFilename, Filename);
  
% return with message indicating invalid function
else
   beep
   disp(['Invalid function name: ' MFileName]);
end
