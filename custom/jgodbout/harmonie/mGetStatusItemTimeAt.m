%MGETSTATUSITEMTIMEAT.M
%
% This function retrieves the status items start time
% given the group and event name at one specific index. 
%
% Usage:
%   [StatusItemTime] = mGetStatusItemTimeAt(Filename, ItemNum, GroupName, EventName)
%
% Input:
%   Filename-       name of the SIG file. Default extension is '.sig'. 
%                   This file must be in MATLAB search path or current working 
%                   directory.
%   ItemNum -       index of the ststus item to be retrieved (1,2, or 3,...)
%   GroupName-      event group name
%   EventName-      event name. If this argument is not specified, any event in 
%                   this groub will count.
%
% Output:
%   StatusItemTime- character string of the start time.
%             
% See also: mGetNumStatusItemsOfEvt, mGetStatusItems, mGetStatusItemTime
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [StatusItemTime] = mGetStatusItemTimeAt(Filename, ItemNum, GroupName, EventName)

StatusItemTime = [];
MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mgetstatusitemtimeat',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs <3 )
       beep
       disp('mGetStatusItemTime requires at least three input arguments.');
       disp('USAGE: [StatusItemTime] = mGetStatusItemTimeAt(Filename, ItemNum, GroupName, EventName)');
       return;
   end
   if (nlhs > 1) 
       beep
       disp('mGetStatusItemTime requires one output argument.');
       disp('USAGE: [StatusItemTime] = mGetStatusItemTimeAt(Filename, ItemNum, GroupName, EventName)');
       return;
   end
   if nrhs == 3,
       EventName = GroupName;
   end
   
   
   % Check to see if the input argument is a string
   if (~ischar(GroupName))
       beep
       disp('mGetStatusItemTimeAt requires the EventName argument to be a character string');
       return;
   end
   if (~ischar(Filename)| size(Filename,1)~=1 )
       beep
       disp('mGetStatusItemTimeAt requires that Sig Filename be a character string.');
       return;
   end
   
   % check if the ItemNum is numeric
   if (~isnumeric(ItemNum))
       beep;
       disp('The input argument ItemNum must be a number!');
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
   [StatusItemTime] =  mSigFileInterFace(MFilename, Filename, ItemNum, GroupName, EventName);
   
% return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
