%MGETSTATUSITEMTIME.M
%
% This function retrieves all the status items start time
% given the group and event name. If no event name is given, all the event 
% item in that group will be retrieved.
%
% Usage:
%       [StatusItemTime] = mGetStatusItemTime(Filename, GroupName, EventName)
%
% Input:
%   Filename-       name of the SIG file. Default extension is '.sig'. 
%                   This file must be in MATLAB search path or current working 
%                   directory.
%   GroupName-      event group name
%   EventName-      event name
%
% Output:
%   StatusItemTime- [numItems x 1] cell array of item start time.
%
%             
% See also: mGetNumStatusItemsOfEvt, mGetStatusItems, mGetStatusItemTimeAt
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [StatusItemTime] = mGetStatusItemTime(Filename, GroupName, EventName)

StatusItemTime = [];
MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mgetstatusitemtime',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs <2 )
       beep
       disp('mGetStatusItemTime requires at least two input arguments.');
       disp('USAGE: [StatusItemTime] = mGetStatusItemTime(Filename, GroupName, EventName )');
       return;
   end
   if (nlhs > 1) 
       beep
       disp('mGetStatusItemTime requires one output argument.');
       disp('USAGE: [StatusItemTime] = mGetStatusItemTime(Filename, GroupName, EventName)');
       return;
   end
   if nrhs == 2,
       EventName = GroupName;
   end
   
   
   % Check to see if the input argument is a string
   if (~ischar(GroupName))
       beep
       disp('mGetStatusItemTime requires the GroupName argument to be a character string');
       return;
   end
   if (~ischar(EventName))
       beep
       disp('mGetStatusItemTime requires the EventName argument to be a character string');
       return;
   end   
   if (~ischar(Filename)| size(Filename,1)~=1 )
       beep
       disp('mGetStatusItemTime requires that Sig Filename be a character string.');
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
   [StatusItemTime] =  mSigFileInterFace(MFilename, Filename, GroupName, EventName);
   
   % return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
