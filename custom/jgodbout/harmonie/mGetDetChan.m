%MGETDETCHAN.M
%
% Given the event (group) name, this function returns the 
% channel names on which the events were detected
%
% Usage:
%      [DetChan] = mGetDetChan(Filename, EvtGrpName)
% 
% Input:
%       Filename-   name of the SIG file. Default extension is '.sig'. 
%                   This file must be in MATLAB search path or current 
%                   working directory.
%       EvtGrpName- event group name
% Output:
%       DetChan-    names of channels on which the events were detected
%                   (cell array)
%  
%
% See also: mGetMtgList, mGetMtgChanFormula
%
% Last modified: Sep. 01, 2004
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DetChan] =mGetDetChan(Filename, EvtName)

Mtgname=[];
MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mgetdetchan',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 2)
       beep
       disp('mGetDetChan requires two input arguments.');
       disp('USAGE: [DetChan] = mGetDetChan(Filename, EventName )');
       return;
   end
   if (nlhs > 1)
       beep
       disp('mGetDetChan requires two input arguments.');
       disp('USAGE: [DetChan] = mGetDetChan(Filename, EventName )');
       return;
   end
   
   % Check to see if the input argument is a string
   if (~ischar(Filename)| size(Filename,1)~=1 )
       beep
       disp('mGetDetChan the Filename argument to be one character string');
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
   
   if (~ischar(EvtName) )
       beep
       disp('mGetDetChan requires EventName argument to be a character string.');
       return;
   end

   
   % Call the MEX function
   [DetChan] = mSigFileInterFace(MFilename, Filename, EvtName);
   
   % return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
