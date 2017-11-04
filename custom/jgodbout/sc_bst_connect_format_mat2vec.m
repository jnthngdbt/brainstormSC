function TFv = sc_bst_connect_format_mat2vec(TFm)
% Converts connectivity data TFm from a matrix form to a compressed  
% vectorized form (upper triangular). Input TFm is of size NcxNcxNf. Output
% TFv is of size Nyx1xNf, where Ny is the number of combinations defined
% from the number of channels Nc: Ny = Nc.(Nc+1)/2.
% 
%       TFv = sc_bst_connect_format_mat2vec(TFm)
% 
% The output TFv is of compressed vectorized form. For expanded vetorized
% form, simply use 
% 
%       TFv = reshape(TFm(:), size(TFm,1)^2, 1, [])
% 
% Example:
% 
%     >> TFm = repmat([1,2;2,3],[1,1,2])
%     TFm(:,:,1) =
%          1     2
%          2     3
%     TFm(:,:,2) =
%          1     2
%          2     3
%     >> TFv = sc_bst_connect_format_mat2vec(TFm)
%     TFv(:,:,1) =
%          1
%          2
%          3
%     TFv(:,:,2) =
%          1
%          2
%          3
% 
% 
% Jonathan Godbout, 2013

[Nc,Nc,Nf] = size(TFm);
Nx = Nc*(Nc+1)/2;

TFv = zeros(Nx,1,Nf);
for i=1:Nf
    TFmi = TFm(:,:,i);
    TFv(:,:,i) = TFmi(triu(ones(Nc))>0);
end

end