function varargout = sc_class_marker( varargin )
% Pseudo-class file for MARKER class. A marker MRK is a simple struct whose
% fields corresponds to features. The only condition is that ALL fields
% (features) must contain the SAME NUMBER OF ELEMENTS IN DIM==2. 
% 
% For example, if MRK contains sleep spindles markers, it could be
% something like this:
% 
%   MRK.amplitude  = [12.4 15.3 11.3 16.3]
%   MRK.position_i = [1754 2585 3574 4725] % Position in samples
%   MRK.position   = [6.8516 10.0977 13.9609 18.4570] % In sec (FS=256Hz)
%   MRK.topography = [10x4 double] % One value for each of 10 channels
%
% Author: Jonathan Godbout, 2013

macro_methodcall;
end

% =========================================================================
% MARKER CLASS METHODS
% =========================================================================

% -------------------------------------------------------------------------
function N = Numel(mrk)

fn = fieldnames(mrk);
if isempty(fn), N = 0; return; end
v = mrk.(fn{1});
if isempty(v)
    N = 0;
else
    if size(v,2)==1
        N = numel(v);
    else
        N = size(v,2);
    end
end

end

% -------------------------------------------------------------------------
function mrk = Add(mrk,mrki)
% If MRK is empty, features correspond to the features of the
% added entity. Otherwise, features considered are those already present in
% MRK; other features of the input entity will be ignored.

if isempty(mrk) || Numel(mrk)==0, mrk = mrki; return; end
   
names = intersect(fieldnames(mrk),fieldnames(mrki));
names_to_remove = setdiff(fieldnames(mrk),fieldnames(mrki));
    
for ii = 1:numel(names)

    if ~isfield(mrki,names{ii}), continue; end

    v = [mrki.(names{ii})]; % 2nd dimension must be observations
    if ischar(v), v = {v}; end
    mrk.(names{ii}) = [mrk.(names{ii}) v];
end  

for ii = 1:numel(names_to_remove)        
    mrk = rmfield(mrk,names_to_remove{ii});
end

end

% -------------------------------------------------------------------------
function [mrko, idx] = Remove(mrki,varargin)
% Remove subset of markers from indexes (logical or not)

[mrko,idx] = Get(mrki,varargin{:});
[mrko,idx] = Get(mrki,~idx);

end

% -------------------------------------------------------------------------
function mrko = Flag(mrki,ref,extend,name)
% Flag markers that concord (position overlap) with another marker set REF.
% The resulting number corresponds to the amount of REF markers that
% concord.
% 
% Both marker sets must have fields POSITION_I, START_I and END_I

if nargin<3, extend = 'both'; end
if nargin<4, name = 'flag_i'; end

mrko = mrki;
switch lower(extend)
    case 'both'
        % Keep both sets with their length
    case 'second'
        mrki.start_i = mrki.position_i;
        mrki.end_i   = mrki.position_i;
    case 'first'
        ref.start_i = ref.position_i;
        ref.end_i   = ref.position_i;
    otherwise
        error('Invalid EXTEND parameter')
end

N = Numel(mrki);
flag_i = zeros(1,N);
for i=1:N
    flag_i(i) = sum((mrki.start_i(i) <= ref.end_i) & (mrki.end_i(i) >= ref.start_i));
end
mrko.(name) = flag_i;

end

% -------------------------------------------------------------------------
function [mrko,idx] = Get(mrki,varargin)
% Get subset of markers from indexes (logical or not)

mrko = struct;
if numel(varargin)==1
    if isnumeric(varargin{1}) || islogical(varargin{1})
        idx = varargin{1};
    end
elseif numel(varargin)==2
    name = varargin{1};
    value = varargin{2};
    idx = strcmpi(value,mrki.(name));
end

names  = fieldnames(mrki);
for ii=1:numel(names)
    ft = mrki.(names{ii});
    if size(ft,2)==1 %%%%%%
        mrko.(names{ii}) = ft; %%%%%%%
    elseif numel(size(ft))==2
        mrko.(names{ii}) = ft(:,idx);
    elseif numel(size(ft))==3
        mrko.(names{ii}) = ft(:,idx,:);
    elseif numel(size(ft))==4
        mrko.(names{ii}) = ft(:,idx,:,:);
    end
end

end

function [m,fidx] = Matrix(mrk,n)
% Return matrix M of features (fields) of marker MRK. Feature names are
% specified by N, which can be cell array of strings of indexes.

if ~exist('n','var')
    n = [];
end

if ~isempty(n)
    if iscell(n)
        Nf = numel(n);
        fidx = zeros(1,Nf);
        fnames = fieldnames(mrk);
        for ii=1:Nf
            fidx(ii) = find(strcmpi(fnames,n{ii}));
        end
    elseif islogical(n)
        fidx = find(n);
        Nf = numel(fidx);
    elseif isnumeric(n)
        fidx = n;
        Nf = numel(fidx);
    else
        error('?');
    end
else
    fnames = fieldnames(mrk);
    fidx = [];
    for ii=1:numel(fnames)
        fi = mrk.(fnames{ii});
        if isnumeric(fi) && size(fi,1)==1
            fidx = [fidx ii];
        end
    end
    Nf = numel(fidx);
end

m = zeros(Nf,Numel(mrk));
for ii=1:Nf
    m(ii,:) = mrk.(fnames{fidx(ii)});
end


end

function [mrk,sort_i] = Sort(mrk,field,mode)
% Reorder the elements of all features (fields) of MRK with the order of
% ordered field named FIELD. Optional is sort mode ('ascend' or 'descend').

if nargin<3, mode = 'ascend'; end

sort_i = 1:Numel(mrk);
if ~isfield(mrk,field)
    warning(['Could not sort on field [' field '] since not available']);
    return
end

v = mrk.(field);
if size(v,1)>1, error('Sorting field must be a 1D feature'); end;
[v_sort, sort_i] = sort(v(:),1,mode);

% Sort features
names = fieldnames(mrk);
for ii=1:numel(names)
    f = mrk.(names{ii});
    if size(f,2)>1
        mrk.(names{ii}) = f(:,sort_i,:);
    end
end

end

