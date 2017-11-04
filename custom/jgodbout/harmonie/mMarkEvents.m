%MMARKEVENTS.M
%
% This function marks events to sig file. It creates a status
% group if it does not exist, otherwise adds new events to an already
% existing group. The event can be instantanesous or interval. Adding to
% an already existing group requires event be compatible, i.e. an 
% instantanesous event can NOT be added to a group of interval events, 
% and vice versa. 
%
% To view the events in Harmonie Reviewer, use mFileClose to close the SIG
% file and save the events before opening it in Reviewer.
% 
% Usage:
%   mMarkEvents(Filename, GrpName, EvtName, MtgName, Chan, EvtColor, ...
%       EvtStartSamp, FeatureName, FeatureValue, EvtDur)
%
% Inputs:
%   Filename-    file to be processed. Default extension is '.sig'. 
%                Must be a character string.
%   GrpName-     event group name. Must be a character string.
%   EvtName-     event name. Must be a character string.
%   MtgName-     montage name. Must be a character string. It is case sensitive.
%                TO MARK TO ALL MONTAGES, input the empty character string '' 
%                or 'ALL'. All the cahnnels will be marked regardless which 
%                channel specified by 'Chan'.
%   Chan-        channel index number. Must be an integer. The channel index 
%                is 1 based. 
%                If this argument is zero, all the channels will be marked.
%   EvtColor-    a character string specifying the color of event marker (not 
%                case sensitive).
%                The supported colors are: MAGENTA, GREEN, CYAN, RED, BLUE, 
%                YELLOW, ORANGE, GRAY
%                LTGRAY, LTBLUE, WHITE, BLACK. Other color names will be taken 
%                as BLACK.
%   EvtStartSamp-vector of start samples where the event will be marked. 
%   FeatureName- a character string containing the property names of the 
%                event. Each property name must be seperated (followed )
%                by a dollar sign '$' in this string. Ex. 'width$amplitude$'.
%                Note: do not put '$' as the first character.
%   FeatureValue-[numEvent x numProperty] numeric matrix corresponding to each 
%                property in FeatureName (in column) and each event (in row).
%   EvtDur-      vector of event durations in samples corresponding to 
%                the start samples. If only a single duration is specifies, 
%                it is assumed the duration is the same for each event to 
%                be marked. Default value is ZERO (instanteneous event). 
%                Note: instanteneous event and interval event can not be 
%                      mixed in one event group.
%
%
% Last modified: Sep 14, 2004
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function mMarkEvents(Filename, GrpName, EvtName, MtgName, Chan, ...
         EvtColor, EvtStartSamp, FeatureName, FeatureValue, EvtDur)


MFilename = lower(mfilename);
EvtColor = upper(EvtColor);

% check the number of the arguments
if (nargin >10 | nargin <9)
   % invalid number of arguments
   beep
   disp('Invalid number of arguments');
   return;
end
if (nargin == 9)
    EvtDur = 0;
end

% Check to see if the input argument is a string
if (~ischar(Filename)| size(Filename,1)~=1 )
      beep
      disp('mGetDataTime requires that Sig Filename be a character string.');
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
    disp('mMarkEvents requires GrpName argument to be a character string.');
    return;
end
if (~ischar(EvtName) )
    beep
    disp('mMarkEvents requires EvtName argument to be a character string.');
    return;
end
if (~ischar(MtgName) )
    beep
    disp('mMarkEvents requires MtgName argument to be a character string.');
    return;
end
if(isempty(MtgName) | strcmp(upper(MtgName),'ALL'))
    MtgName='ALL'; % all montages
    Chan =0 ;      % all channels
end
if (~ischar(EvtColor) )
    beep
    disp('mMarkEvents requires EvtColor argument to be a character string.');
    return;
end
if (~ischar(FeatureName)|isempty(FeatureName) )
    beep
    disp('mMarkEvents requires FeatureName argument to be a non-empty character string.');
    return;
end

% check numeric arguments
if (~isnumeric(Chan))
    beep
    disp('mMarkEvents requires Chan argument to be a number.');
    return;
end
if (~isnumeric(EvtStartSamp))
    beep
    disp('mMarkEvents requires EvtStartSamp argument be numberic.');
    return;
end
if (~isnumeric(EvtDur))
    beep
    disp('mMarkEvents requires EvtDur argument be numberic.');
    return;
end

% change the numeric range
if (~strcmp(MtgName,'ALL'))
    ChanFormula = mGetMtgChanFormula(Filename, MtgName);
    if(isempty(ChanFormula))
        beep
        disp(['The montage ' MtgName ' does NOT exist!']);
        return;
    end
    if (Chan<0 | Chan>size(ChanFormula,1))
        beep
        disp('Channel number is out of range!');
        return;
    end
end

[m n]=size(EvtStartSamp);
if m<n
    EvtStartSamp=EvtStartSamp';
end
[m n]=size(EvtDur);
if m<n
    EvtDur=EvtDur';
end
if(size(EvtDur,2)>1 | size(EvtStartSamp,2)>1)
    beep
    disp('mMarkEvents requires EvtDur and EvtStartSamp arguments be 1 dimensional number array.');
    return;
end
if (length(EvtDur) > length(EvtStartSamp))
    EvtDur = EvtDur(1:length(EvtStartSamp));
end
if(length(EvtDur) < length(EvtStartSamp))
    EvtDur=EvtDur(1)*ones(1,length(EvtStartSamp));
end
if (~isempty(find(EvtStartSamp<0)))
    beep
    disp('EvtStartSamp can not be negative!')
    return;
end
if (~isempty(find(EvtDur<0)))
    beep
    disp('EvtDur can not be negative!')
    return;
end        

% check the duration compatibility
if (min(EvtDur)==0 & max(EvtDur >0))
    beep
    disp('Instanteneous event and interval event can not be mixed. Please check EvtDur argument.');
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
    
if length(EvtStartSamp)~=size(FeatureValue,1)
    beep
    disp('Number of FeatureValue does NOT match number of events (EvtStartSamp).')
    return;
end

if length(findstr(FeatureName, '$'))~=size(FeatureValue,2)
    beep
    disp('Number of FeatureValue does NOT match number of properties in FeatureName.')
    return;
end

mSigFileInterFace(MFilename, Filename, GrpName, EvtName, MtgName, Chan, ...
         EvtColor, EvtStartSamp, FeatureName, FeatureValue, EvtDur);




