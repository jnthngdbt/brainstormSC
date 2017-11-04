%MGETMTGCHANFORMULA.M
%
% This function extracts the reformatting montage 
% formula given the filename and the montage name.
%
% Usage:
%   [ChanFormula]= mGetMtgChanFormula(Filename, MtgName)
%
% Input:
%   Filename-       name of the SIG file. Default extension is '.sig'. 
%                   This file must be in MATLAB search path or current working 
%                   directory.
%   MtgName -       montage name. It is case sensitive.
%
% Output:
%   ChanFormula-   [NumChannel x 1] cell array of the channel formula
%             
%
% See also: mGetMtgList, mGetDetMtg, mGetChanLabel
%
% Last modified: Dec. 18, 2002
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ChanFormula] =mGetMtgChanFormula(Filename, MtgName)

ChanFormula = [];
MFilename = lower(mfilename);
%disp(MFilename);

if MFilename == 'mgetmtgchanformula',
   % Check if correct number of input/output arguments
   nrhs = nargin;
   nlhs = nargout;
   if (nrhs ~= 2)
       beep
       disp('mGetMtgChanFormula requires two input arguments.');
       disp('USAGE: [ChanFormula] = mGetMtgChanFormula(Filename, MtgName )');
       return;
   end
   if (nlhs > 1)
       beep
       disp('mGetMtgChanFormula requires one output arguments.');
       disp('USAGE: [ChanFormula] = mGetMtgChanFormula(Filename, MtgName )');
       return;
   end
   
   % Check to see if the input argument is a string
   if (~ischar(Filename)| size(Filename,1)~=1 )
       beep
       disp('mGetMtgChanFormula the Filename argument be a character string');
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
   
   if (~ischar(MtgName) )
       beep
       disp('mGetMtgChanFormula requires MtgName argument be a character string.');
       return;
   end
   
   % Call the MEX function
   [ChanFormula] = mSigFileInterFace(MFilename, Filename, MtgName);
   
%    % Reformat the ChanFormula data
%    if isempty(ChanFormulaTmp)
%        beep
%        ChanFormula=[];
%        disp(['Montage (' MtgName ') does NOT exist!']);
%        return;
%    end
%    
%     ChanFormulaTmp = [ChanFormulaTmp '$'];
%     Pos = findstr(ChanFormulaTmp, '$');
%     N=max(diff(Pos));
%     Nchan = length(Pos);
%     for ij =1:Nchan-1
%        x = ChanFormulaTmp(Pos(ij)+1:Pos(ij+1)-1);
%        x(N) = ' ';
% 	   ChanFormula =[ChanFormula; x];
%     end
    
   % return with message indicating invalid function
else
    beep
    disp(['Invalid function name: ' MFileName]);
end
