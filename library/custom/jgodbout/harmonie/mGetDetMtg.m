%MGETDETMTG.M
%
% Given the event (group) name, this function returns the 
% corresponding montage name.
%
% Usage:
%      [MtgName] = mGetDetMtg(Filename, EvtGrpName)
% 
% Input:
%       Filename-   name of the SIG file. Default extension is '.sig'. 
%                   This file must be in MATLAB search path or current 
%                   working directory.
%       EvtGrpName- event group name
% Output:
%       MtgName-    Montage name associated with EvtName. If the event 
%                   item is associated to all the monmtages, the rerurn 
%                   will be 'All Mtgs'.
%
% Note: For 6.0 version, the montage from this function is the first
%       montage in the montage list that contains all the channles on which
%       detections are marked. It may not be the montage that was used to
%       detect the events.
%
% See also: mGetMtgList, mGetMtgChanFormula
%
% Last modified: Apr. 21, 2005
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [MtgName] =mGetDetMtg(Filename, EvtName)

Mtgname=[];
MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mgetdetmtg',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 2)
       beep
       disp('mGetDetMtg requires two input arguments.');
       disp('USAGE: [MtgName] = mGetDetMtg(Filename, EventName )');
       return;
   end
   if (nlhs > 1)
       beep
       disp('mGetDetMtg requires one output argument.');
       disp('USAGE: [MtgName] = mGetDetMtg(Filename, EventName )');
       return;
   end
   
   % Check to see if the input argument is a string
   if (~ischar(Filename)| size(Filename,1)~=1 )
       beep
       disp('mGetDetMtg the Filename argument to be one character string');
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
       disp('mGetDetMtg requires EventName argument to be a character string.');
       return;
   end

   
   [DetChan] = mGetDetChan(Filename, EvtName);
   
   if size(DetChan,1) == 0
       MtgName = 'All Mtgs';
       return;
   end
       
   
   [MtgList] = mGetMtgList(Filename);
   
   for i=1:size(MtgList,1)
       ChanLabel = mGetChanLabel(Filename, char(MtgList(i)));
       bEx = 1;
       for j=1:size(DetChan,1)
           bFound = 0;
           for k=1:size(ChanLabel,1)
               if strcmp(char(ChanLabel(k)), char(DetChan(j)))
                   bFound = 1;
                   break;
               end
           end
           if bFound == 0
               bEx = 0;
               break;
           end
       end
       
       if bEx == 1
           break;
       end
   end
   if (bEx == 1)
       MtgName = char(MtgList(i));
   else
       MtgName = 'All Mtgs';
   end
       
                   
   
   % Call the MEX function
   % [MtgName] = mSigFileInterFace(MFilename, Filename, EvtName);
   
   % return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
