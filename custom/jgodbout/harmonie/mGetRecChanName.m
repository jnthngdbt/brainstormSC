%MGETRECCHANNAME.M
%
% Given the file name, this function returns the electrode 
% names of the recording montage.
%
% Usage:
%       [RecChanNames] = mGetRecChanName(Filename)
%
% Input:
%   Filename-     name of the SIG file. Default extension is '.sig'. 
%                 This file must be in MATLAB search path or current working 
%                 directory.
% Output:
%   RecChanNames- [numChannel x 1] cell array of channel names. 
%             
%
% See also: mGetNumChan, mGetMtgList, mGetChanLabel, mGetMtgChanFormula
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [RecChanNames] =mGetRecChanName(Filename)

RecChanNames = [];
MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mgetrecchanname',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 1)
       beep
       disp('mGetRecChanName requires one input arguments.');
       disp('USAGE: [RecChanNames] = mGetRecChanName(Filename)');
       return;
   end
   if (nlhs > 1)
       beep
       disp('mGetRecChanName requires one output arguments.');
       disp('USAGE: [RecChanNames] = mGetRecChanName(Filename)');
       return;
   end
   
   % Check to see if the input argument is a string
   if (~ischar(Filename)| size(Filename,1)~=1)
       beep
       disp('mGetRecChanName requires the Filename argument be a character string');
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
    [RecChanNames] = mSigFileInterFace(MFilename, Filename);
    

    % return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
