%MCHANGESTATUSITEMEVTNAME.M
%
% This function changes the event name within one event group. It is used to calssify 
% one event group into sub-groups.
%
% Usage:
%      mChangeStatusItemEvtName(Filename, EvtGrpName, EvtChan, EvtBeginSamp, EvtName)
%
% Input:
%      Filename-       SIG File including the path (default extension: .SIG)
%      EvtGrpName-     Name of the event group
%      EvtChan-        Event Channel Name (string)
%      EvtBegingSamp-  Beginning sample number of the event 
%      EvtName-        Name to be assigned to the event
%
% Here, together with event group name, the event channel number and sample number 
% define the event item which will be assigned the new name (EvtName). It can  
% apply to single event. For multiple events, use it repeatedly.
%
% see also: mGetStatusItems, mGetNumStatusItemOfEvt
%
% Last modified: Sep. 01, 2004
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mChangeStatusItemEvtName(Filename, EvtGrpName, EvtChan, EvtBeginSamp, EvtName)

MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mchangestatusitemevtname',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 5)
      disp('mGetStatusItems requires five input arguments.');
      disp('USAGE: mChangeStatusItemEvtName(Filename, EvtGrpName, EvtChan, EventBeginSamp, EventName)');
      return;
   end
   
   % Check to see if the input argument is a string
   if (~ischar(Filename)| size(Filename,1)~=1 )
        beep
        disp('mChangeStatusItemEvtName the Filename argument to be one character string');
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
%    Filename = which(Filename);
   
   % check the file R/W status
   fclose all;
%    fidG=fopen(Filename,'a');
%    fidS=fopen([P '\' N '.sts'],'a');
%    fclose all;
%    if(fidG==-1)|(fidS==-1)
%        beep
%        disp(['File ' Filename ' can NOT be opened for writing!']); 
%        return;
%    end
   

   if (~ischar(EvtGrpName))
       beep
       disp('mChangeStatusItemEvtName requires the EventGroupName argument to be a character string');
       return;
   end
   if (~ischar(EvtName) )
       beep
       disp('mChangeStatusItemEvtName requires EventName argument to be a character string.');
       return;
   end
   if (~ischar(EvtChan) )
       beep
       disp('mChangeStatusItemEvtName requires EventChan argument to be a character string.');
       return;
   end
   
   if (~isnumeric(EvtBeginSamp) )
       beep
       disp('mChangeStatusItemEvtName requires EventBeginSamp argument to be numeric.');
       return;
   end
   
   if (length(EvtChan)~=length(EvtBeginSamp))
       beep
       disp('EvtChan and EvtBeginSamp must be of same size!');
       return;
   end

   
   % Call the MEX function 
   mSigFileInterFace(MFilename, Filename, EvtGrpName, EvtName, EvtChan, EvtBeginSamp);
   
   % return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
