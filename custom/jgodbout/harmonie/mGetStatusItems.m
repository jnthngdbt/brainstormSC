%MGETSTATUSITEMS.M
%
% This function retrieves all the status items and their properties 
% given the group name and event name. The properties names are described
% by the PropHeader. If no event name is given, all the event item in that
% group will be retrieved.
%
% Usage:
%       [PropHeader, StatusItems, DetChan] = mGetStatusItems(Filename, GroupName, EventName)
%
% Input:
%   Filename-   name of the SIG file. Default extension is '.sig'. 
%               This file must be in MATLAB search path or current working 
%               directory.
%   GroupName-  event group name
%   EventName-  event name
%
% Output:
%   PropHeader- property header decribing properties corresponding to the 
%               columns in StatusItems.
%   StatusItems-property values corresponding to the properties given in
%               PropHeader. 
%               Note: all the index numbers follows C convention, i.e. 0 based.
%   DetChan-    names of channels on which the events were detected
%                   (cell array)
%   
%
% See also: mGetNumStatusItemsOfEvt, mGetStatusItemTime, mGetStatusItemTimeAt
%
% Last modified: Apr. 20, 2005
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PropHeader, StatusItems, DetChan] = mGetStatusItems(Filename, GroupName, EventName)

PropHeader = [];
StatusItems = [];
MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mgetstatusitems',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs <2 )
       beep
       disp('mGetStatusItems requires at least two input arguments.');
       disp('USAGE: [PropHeader, StatusItems, DetChan] = mGetStatusItems(Filename, GroupName, EventName )');
       return;
   end
   if (nlhs > 3) 
       beep
       disp('mGetStatusItems requires three output argument.');
       disp('USAGE: [PropHeader, StatusItems, DetChan] = mGetStatusItems(Filename, EventName )');
       return;
   end
   if nrhs == 2,
       EventName = 'AllEventInGroup';
   end
   
   
   % Check to see if the input argument is a string
   if (~ischar(GroupName)| size(GroupName,1)~=1)
       beep
       disp('mGetStatusItems requires the EventName argument to be a character string');
       return;
   end
   
   if (~ischar(Filename) | size(Filename,1)~=1)
       beep
       disp('mGetStatusItems requires that Sig Filename be a character string.');
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
  [PropHeader, StatusItems, DetChan] =  mSigFileInterFace(MFilename, Filename, GroupName, EventName);
   
% return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
