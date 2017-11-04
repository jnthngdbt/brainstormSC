%MGETEVTGRPLIST.M
%
% Given the file name, this function returns the event group list.
%
% Usage:
%   [GrpList] = mGetEvtGrpList(Filename)
%
% Input:
%   Filename-  name of the SIG file. Default extension is '.sig'. 
%              This file must be in MATLAB search path or current working 
%              directory.
% Output:
%   GrpList -  [numGrps x 1] cell array of group names. 
%             
%
% See also: mGetStatusItems.m 
%
% Created: Nov. 17, 2003
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [GrpList] = mGetEvtGrpList(Filename)


MFilename = lower(mfilename);
%disp(MFilename);


if MFilename == 'mgetevtgrplist',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 1)
      beep
      disp('mGetEvtGrpList requires one input arguments.');
      return;
   end
	if (nlhs > 1) 
      beep
      disp('mGetEvtGrpList requires 1 output argument.');
      return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(Filename) | size(Filename,1)~=1 )
     beep
     disp('mGetEvtGrpList requires that Sig Filename be one character string.');
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
  [GrpList] =  mSigFileInterFace(MFilename, Filename);
  
% return with message indicating invalid function
else
   beep
   disp(['Invalid function name: ' MFileName]);
end
