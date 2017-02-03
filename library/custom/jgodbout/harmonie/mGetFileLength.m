%MGETFILELENGTH.M
%
% This function returns the size of a SIG file in 
% number of records, samples, and seconds.
%
% Usage:
%   [NumRecs, NumSamps, NumSecs] = mGetFileLength(Filename)
%
% Input:
%   Filename-  name of the SIG file. Default extension is '.sig'. 
%              This file must be in MATLAB search path or current working 
%              directory.
% Output:
%   NumRecs -  number of records
%   NumSamps-  number of samples
%   NumSecs -  number of seconds
%             
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [NumRecs, NumSamps, NumSecs] = mGetFileLength(Filename)

NumOfChan = 0;
MFilename = lower(mfilename);
%disp(MFilename);


if MFilename == 'mgetfilelength',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 1)
      beep
      disp('mGetFileLength requires one input arguments.');
      return;
   end
	if (nlhs > 3) 
      beep
      disp('mGetFileLength requires 3 output argument.');
      return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(Filename) | size(Filename,1)~=1 )
     beep
     disp('mGetFileLength requires that Sig Filename be one character string.');
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
  [FileLength] =  mSigFileInterFace(MFilename, Filename);
  
  NumRecs = fix(FileLength(1));
  NumSamps = fix(FileLength(2));
  NumSecs = FileLength(3);
  
% return with message indicating invalid function
else
   beep
   disp(['Invalid function name: ' MFileName]);
end
