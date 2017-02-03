%MGETRECSTARTTIME.M
%
% This function extract the start time (i.e. the 
% time of the first recorded sample in the given 
% SIG file. Output is a string.
%
% Usage:
%   [RecordingStartTime] = mGetRecStartTime(Filename)
%
% Input:
%   Filename-           name of the SIG file. Default extension is '.sig'. 
%                       This file must be in MATLAB search path or current 
%                       working directory.
% Output:
%   RecordingStartTime- start time of the recording
%             
%
% See also: mGetDataTime
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RecordingStartTime] = mGetRecStartTime(Filename)

RecordingStartTime = ' ';
MFilename = lower(mfilename);
%disp(MFilename);


if MFilename == 'mgetrecstarttime',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 1)
       beep
       disp('mGetRecStartTime requires one input arguments.');
       return;
   end
   if (nlhs > 1) 
       beep
       disp('mGetRecStartTime requires one output argument.');
       return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(Filename)| size(Filename,1)~=1 )
      beep
      disp('mGetRecStartTime requires that Sig Filename be a character string.');
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
  [RecordingStartTime] =  mSigFileInterFace(MFilename, Filename);
  
% return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
