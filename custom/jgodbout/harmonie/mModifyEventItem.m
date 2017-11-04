%MMODIFYEVENTITEM.M
%
% This function modififies propertoes of a specific event item to sig file. 
% The event can be instantanesous or interval.  If event name is empty, all
% the event in this group will be counted.
%
% To view the events in Harmonie Reviewer, use mFileClose to close the SIG
% file and save the events before opening it in Reviewer.
% 
% Usage:
%   mModifyEventItem(Filename, GrpName, EvtName, EvtIdx, FeatureName, FeatureValue)
%
% Inputs:
%   Filename-    file to be processed. Default extension is '.sig'. 
%                Must be a character string.
%   GrpName-     event group name. Must be a character string.
%   EvtName-     event name. Must be a character string.
%   EvtIdx-      the event index to be modified given the Grpname and EvtName (0 based)
%   FeatureName- a character string containing the property names of the 
%                event. Each property name must be seperated (followed )
%                by a dollar sign '$' in this string. Ex. 'width$amplitude$'.
%                Note: do not put '$' as the first character.
%   FeatureValue-[1 x numProperty] numeric row vector corresponding to each 
%                property in FeatureName.
%
% See also: mGetStatusItems, mGetStatusItemTimeAt
%
% Last modified: Nov 10, 2005
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mModifyEventItem(Filename, GrpName, EvtName, EvtIdx, FeatureName, FeatureValue)


MFilename = lower(mfilename);

% check the number of the arguments
if (nargin ~= 6)
   % invalid number of arguments
   beep
   disp('Invalid number of arguments');
   return;
end

% Check to see if the input argument is a string
if (~ischar(Filename)| size(Filename,1)~=1 )
      beep
      disp('mModifyEventItem requires that Sig Filename be a character string.');
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

% check the file R/W status
fclose all;
% fidG=fopen(Filename,'a');
% fidS=fopen([P '\' N '.sts'],'a');
% fclose all;
% if(fidG==-1)|(fidS==-1)
%     beep
%     disp(['File ' Filename ' can NOT be opened for writing!']); 
%     return;
% end

% check character string arguments
if (~ischar(GrpName) )
    beep
    disp('mModifyEventItem requires GrpName argument to be a character string.');
    return;
end

if isempty(EvtName)
    EvtName = GrpName;
end
if (~ischar(EvtName) )
    beep
    disp('mModifyEventItem requires EvtName argument to be a character string.');
    return;
end

if (~ischar(FeatureName)|isempty(FeatureName) )
    beep
    disp('mModifyEventItem requires FeatureName argument to be a non-empty character string.');
    return;
end

% check numeric arguments
if (~isnumeric(EvtIdx) | EvtIdx<0)
    beep
    disp('mModifyEventItem requires EvtIdx argument to be a positive number.');
    return;
end


if strcmp(FeatureName(1),'$')
   FeatureName(1)=[];
end

if isempty(FeatureName)
    beep
    disp('No property is specified in the argument FeatureName!');
    return;
end
if ~strcmp(FeatureName(end),'$')
   FeatureName(end+1)='$';
end


if length(findstr(FeatureName, '$'))~=length(FeatureValue)
    beep
    disp('Number of FeatureValue does NOT match number of properties in FeatureName.')
    return;
end

mSigFileInterFace(MFilename, Filename, GrpName, EvtName, EvtIdx, FeatureName, FeatureValue);




