function sc_bst_install
% All files (of editable extensions, but mainly M-files) in folders whose 
% root is the folder of the present function will be copied in Brainstorm 
% folder, respecting relative paths. It also respected modification date; 
% newer file overwrites older.
% 
% Jonathan Godbout, 2013

dir0 = cd;

% Extensions of files to link
link_ext = {'.m','.dll','.p','.mexa64','.mexmaci64','.mexw32','.mexw64','.sel'};

% Get current function file name parts
[link_path,link_name,tmp] = fileparts(which(mfilename('fullpath')));
if ~isempty(link_path), link_path = [link_path,filesep];
else                    link_path = [cd,filesep];
end

% Add toolbox function to path, then remove subfolders. The subfolders
% contain functions that are linked (copied) in Brainstorm library. We want
% that Matlab sees only those in Brainstorm folders.
addpath(genpath(link_path));
rmpath(strrep(genpath(link_path),[link_path,pathsep],''));
savepath;

% Scan through subfolders and update files
find_do(link_path,link_ext,@(x)update_bst_function(x,link_path,link_name));

% Make sure that all subdirectories in custom are added (not done by BST)
addpath(genpath([bst_path,'custom']));
savepath;

cd(dir0);

end

function update_bst_function(src,link_path,link_name)
% n is the fullname of the file

if isequal(src,which(link_name)), return; end

% Set the destination (brainstorm) full file name
rel = src;
rel(1:length(link_path)) = [];
dst = [bst_path, rel];

% If folder doesn't exist, create it
[pDst,nDst,eDst] = fileparts(dst);
if ~isdir(pDst), mkdir(pDst); end

% Which of the link or brainstorm file is the newest
dSrc = dir(src);
dDst = dir(dst);
if ~isfield(dSrc,'datenum'), dSrc.datenum = datenum(dSrc.date); end
if ~isfield(dDst,'datenum'), dDst.datenum = datenum(dDst.date); end
if ~isempty(dDst)
    if dDst.datenum==dSrc.datenum
        return;
    elseif dDst.datenum>dSrc.datenum % BST is the newest
        tmp = src; src = dst; dst = tmp; % Switch
    end
end

fprintf('\nUPDATING FILE [%s]\nWITH FILE     [%s]\n',dst,src);
copyfile(src,dst);

end

function p = bst_path
[p,aaa,bbb] = fileparts(which('brainstorm'));
if isempty(p);
    error('Brainstorm not in path!');
end
p = [p,filesep];

end

function find_do(dir0,e,f,varargin)
%% FIND_DO
% Search files of specified type and do something on them
% 
%   FIND_DO(DIR0,E,F)
% 
% DIR0: root folder (ex. cd)
% E: extension of the wanted file (ex.: '.m' or '.png'). If empty, any file
% F: handle to a function (ex.: @open)
% 
% Consider one folder and all its subfolders as a tree, the content of one
% of these folders as leaves and the initial folder as the root. This
% function is an algorithm that starts from the root DIR0 and seeks into
% all its subfolders. In each leaf, it searches for files of extension E.
% If some are present, an operation is applied on them, defined in the
% function whose handle F is sent in the present function.
% 
% In formal word, for every file X of extension E, we call F(X)

init_dir = cd;

cd(dir0);

L = 1;                  % Level
S(L) = 1;               % Step: folder scanned at specific level
F{1} = {[dir0 filesep]};  % Folder name at first level

L = L+1;                % Jump to higher level
S(L) = 1;               % Initialize step
F{L} = {};              % Initialize folder names

while 1                 % We hope that you converge
    nodes = dir;
    
    for ii=1:numel(nodes)
        [pathstr,name,ext] = fileparts(nodes(ii).name);
        if isempty(pathstr), pathstr = [cd filesep]; end
        if nodes(ii).isdir && ~strcmpi(nodes(ii).name,'.') && ~strcmpi(nodes(ii).name,'..')
            F{L}(end+1) = {[F{L-1}{S(L-1)} name ext filesep]};
        else
            if isempty(e) && ~ismember(ext,{'.','..'})
                f([pathstr,name,ext],varargin{:});
            else
                switch lower(ext)
                    case lower(e)
                        f([pathstr,name,ext],varargin{:});
                    otherwise
                end
            end
        end
    end

    if isempty(F{L})
        while 1
            F = F(1:L-1);           % Clear this empty level (no folder)
            S = S(1:L-1);           % No step in an empty level      
            L = L-1;                % Jump back to lower level
            S(L) = S(L)+1;          % Step to next folder
            if L==1, cd(dir0); return; end
            if S(L) <= numel(F{L}), break; end
        end
    end
    cd(F{L}{S(L)})          % Go in the step folder to higher level
    L = L+1;                % Jumped to higher level
    S(L) = 1;               % Initialize step
    F{L} = {};              % Initialize folder names
end

cd(init_dir);
end

% function n = link_name
% n = 'sc_bst_install'; % MAKE SURE IT IS THE GOOD NAME %%%%%%%%%%%%%
% end
% 
% function p = link_path
% [p,~,~] = fileparts(which(link_name)); p = [p,filesep];
% end
