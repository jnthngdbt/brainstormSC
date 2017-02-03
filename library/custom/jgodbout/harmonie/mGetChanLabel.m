%MGETCHANLABEL.M
%
% This function extracts the montage channel
% labels given the filename and the montage name.
%
% Usage:
%   [ChanLabel]= mGetChanLabel(Filename, MtgName)
%
% Input:
%   Filename-   name of the SIG file. Default extension is '.sig'. 
%               This file must be in MATLAB search path or current working 
%               directory.
%   MtgName -   montage name. It is case sensitive.
%
% Output:
%   ChanLabel-  [numChannel x 1] cell array of the channel labels
%             
%
% See also: mGetMtgList, mGetDetMtg, mGetMtgChanFormula, mGetRecChanName
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ChanLabel] =mGetChanLabel(Filename, MtgName)

ChanLabel = [];
MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mgetchanlabel',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 2)
       beep
       disp('mGetChanLabel requires two input arguments.');
       disp('USAGE: [ChanLabel] =mGetChanLabel(Filename, MtgName)');
       return;
   end
   if (nlhs > 1)
       beep
       disp('mGetChanLabel requires one output arguments.');
       disp('USAGE: [ChanLabel] =mGetChanLabel(Filename, MtgName)');
       return;
   end
   
   % Check to see if the input argument is a string
   if (~ischar(Filename)| size(Filename,1)~=1 )
       beep
       disp('mGetChanLabel requires the Filename argument be a character string');
       return;
   end
   % Check the extension (default extension is '.sig')
   [P N E]=fileparts(Filename);
   if isempty(E)
      Filename=[Filename '.sig'];
   end
  
   % check the existence of the file
   if exist(Filename,'file')==0
     beep
     disp(['File ' Filename ' does NOT exist!']); 
     return;
   end   
%    Filename = which(Filename);
   
   if (~ischar(MtgName) )
       beep
       disp('mGetChanLabel requires MtgName argument be a character string.');
       return;
   end
   
   % Call the MEX function
   [ChanLabel] = mSigFileInterFace(MFilename, Filename, MtgName);

   % return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
