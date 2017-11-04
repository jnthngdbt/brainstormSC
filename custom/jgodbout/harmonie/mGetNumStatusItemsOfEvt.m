%MGETNUMSTATUSITEMOFEVT.M
%
% Given the event group name, this function returns the 
% number of event of that type. 
%
% Usage:
%   [NumStatusItems] = mGetNumStatusItemsOfEvt(Filename, EventGrpName)
%
% Input:
%   Filename-       name of the SIG file. Default extension is '.sig'. 
%                   This file must be in MATLAB search path or current 
%                   working directory.
%   EventGrpName-   event group name
%
% Output:
%   NumStatusItems- total number of event items in this group
%
%             
% See also: mGetStatusItems, mGetStatusItemTime, mGetStatusItemTimeAt
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NumStatusItems] = mGetNumStatusItemsOfEvt(Filename, EventGrpName)

NumStatusItems=[];
MFilename = lower(mfilename);

if MFilename == 'mgetnumstatusitemsofevt',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 2)
       beep
       disp('mGetNumStatusItemsOfEvt requires two input arguments.');
       disp('USAGE: [NumStatusItems] = mGetNumStatusItemsOfEvt(Filename, EventGrpName )');
       return;
   end
   if (nlhs > 1) 
       beep
       disp('mGetNumStatusItemsOfEvt requires one output argument.');
       disp('USAGE: [NumStatusItems] = mGetNumStatusItemsOfEvt(Filename, EventGrpName )');
       return;
   end
  
  % Check to see if the input argument is a string
  if (~ischar(EventGrpName))
      beep
      disp('mGetNumStatusItemsOfEvtrequires the EventGrpName argument to be a character string');
      return;
  end
  
  if (~ischar(Filename)| size(Filename,1)~=1 )
     beep
     disp('mGetNumStatusItemsOfEvt requires that Sig Filename be a character string.');
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
  
  % Call the MEX function 
  [NumStatusItems] =  mSigFileInterFace(MFilename, Filename, EventGrpName);
  
% return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
