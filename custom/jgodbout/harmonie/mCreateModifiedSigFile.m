%MCREATEMODIFIEDSIGFILE.M
%
% This function creates a new SIG file by modifying ONE channel
% of an existing SIG file on record base. It supports two types of modifications:
% (1) adding data to the channel; (2) replace the channel with new data.
% The data is from a binary data file, which must be prepared in advance.
% The data format must be 4-byte float (in MATLAB, use fwrite(.,.,'float')).
% This function will overwrite the existing SIG file. The modification 
% begins from the first record after stTime till the end of SIG file 
% or data file, whichever comes first. Fraction of record in data file will 
% be ignored.
% 
% Note: (1) If you want to change multiple channels, you have to prepare multiple
%           data files and run this function multiple times on the same SIG file.
%       (2) This function applies to recording montage.
%
% Usage:
%       mCreateModifiedSigFile(sigFilename, dataFilename, stTime, Chan, Rule)
%
% Inputs:
%       sigFilename- existing SIG file
%       dataFilename-binary data file 
%       stTime-      The start time in the SIG file where modification will
%                    take place. It must be in an array in the form of 
%                    [year, month, day, hour, minute, second]
%       Chan-        the channel to be modified
%       Rule-        type of modification: 
%                       'ADD'-add the data to the existing channel
%                       'REP'-replace the channel with the data
%
% See also: mGetRecStartTime, mGetNumChan, mGetDataTime
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mCreateModifiedSigFile(sigFilename, dataFilename, stTime, Chan, Rule)

MFilename = lower(mfilename);

% check the number of the arguments
if (nargin ~=5)
   % invalid number of arguments
   beep
   disp('Invalid number of arguments');
   return;
end

% Check to see if the input argument is a string
if (~ischar(sigFilename)| size(sigFilename,1)~=1 )
      beep
      disp('mCreateModifiedSigFile requires that sigFilename be a character string.');
      return;
end
  
% Check the extension (default extension is '.sig')
[P N E]=fileparts(sigFilename);
if isempty(E)
      sigFilename=[sigFilename '.sig'];
end
  
% chack the existence of the file
if exist(sigFilename,'file')==0
     beep
     disp(['File ' sigFilename ' does NOT exist!']); 
     return;
end
% Filename = which(Filename);

% check the file R/W status
fclose all;
fidG=fopen(sigFilename,'a');
fidS=fopen([P '\' N '.sts'],'a');
fclose all;
if(fidG==-1)|(fidS==-1)
    beep
    disp(['File ' sigFilename ' can NOT be opened for writing!']); 
    return;
end

% Check to see if the input argument is a string
if (~ischar(dataFilename)| size(dataFilename,1)~=1 )
      beep
      disp('mCreateModifiedSigFile requires that dataFilename be a character string.');
      return;
end
% chack the existence of the file
if exist(dataFilename,'file')==0
     beep
     disp(['File ' dataFilename ' does NOT exist!']); 
     return;
end



% check the stTime argument
if (~isnumeric(stTime) )
      beep
      disp('mCreateModifiedSigFile requires that stTime be a numeric array.');
      return;
end

if length(stTime)~= 6
      beep
      disp('Please specify the time in the form of [year month day hour minute second].');
      return;
end

% check stTime argument
if (~isnumeric(Chan) )
      beep
      disp('mCreateModifiedSigFile requires that Chan be numeric.');
      return;
end

numChan = mGetNumChan(sigFilename);
if (Chan<0 |Chan>numChan-1)
      beep
      disp('The channel number is out of range.');
      return;
end  

% check Rule argument
if (~ischar(Rule))
      beep
      disp('mCreateModifiedSigFile requires that Rule be a character string.');
      return;
end

if (~strcmp(upper(Rule),'ADD'))& (~strcmp(upper(Rule),'REP'))
      beep
      disp('Rule should be either ''ADD'' or ''REP''! ');
      return;
end 

% call MEX function
mSigFileInterFace(MFilename, sigFilename, dataFilename, stTime, Chan, Rule);
  
