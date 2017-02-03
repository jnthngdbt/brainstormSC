function TFm = sc_bst_connect_format_vec2mat(TFv,sym)
% Converts connectivity data TFV from a compressed vectorized form to a
% matrix form. TFV is of size Nyx1xNf, where Ny is the number of 
% combinations defined from the number of channels Nc: Ny = Nc.(Nc+1)/2. 
% Dimension 2 must be of size 1. Output TFm is of size NcxNcxNf.
% 
%       TFm = sc_bst_connect_format_vec2mat(TFv)
%       TFm = sc_bst_connect_format_vec2mat(TFv,sym)
% 
% The optional input 'sym' is the symetry of the matrix ('invert', for 
% negative inverse or 'conjugate', for complex conjugate). By defaults:
% 'conjugate'.
% 
% TFv must be of compressed form (upper triangular). For expanded forms,
% simply use 
% 
%       TFm = reshape(TFv,[Nc,Nc])
% 
% where Nc is the number of channels:
% 
%       Nc = (-1+sqrt(1+ 8*size(TFv,1)))/2
% 
% Example:
% 
%     >> TFv = repmat((1:3)',[1,1,2])
%     TFv(:,:,1) =
%          1
%          2
%          3
%     TFv(:,:,2) =
%          1
%          2
%          3
%     >> TFm = sc_bst_connect_format_vec2mat(TFv)
%     TFm(:,:,1) =
%          1     2
%          2     3
%     TFm(:,:,2) =
%          1     2
%          2     3
% 
% 
% Jonathan Godbout, 2013

if nargin<2, sym = 'conjugate'; end

[Nx,Nt,Nf] = size(TFv);

if Nt~=1, error('Dimension 2 (time) must be of size 1'); end

Nc = (-1+sqrt(1+ 8*size(TFv,1)))/2;

% Extracted from Expand
% Generate all the indices
[iAall,iBall] = meshgrid(1:Nc,1:Nc);
% Find the values below/above the diagonal
[iA,iB] = find(iBall <= iAall);

indAll1 = sub2ind([Nc,Nc], iA(:), iB(:));
indAll2 = sub2ind([Nc,Nc], iB(:), iA(:));

TFm = zeros(Nc,Nc,Nf);
for i=1:Nf
    TFmi = TFm(:,:,i);
    TFmi(indAll1) = TFv(:,:,i);
    switch lower(sym)
        case 'conjugate'
            TFmi(indAll2) = conj(TFv(:,:,i));
        case 'invert'
            TFmi(indAll2) = -TFv(:,:,i);
    end
    TFm(:,:,i) = TFmi;
end

end