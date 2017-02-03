function i = sc_bst_connect_vector_diagind(x,isExpanded)
% Returns the indexes, of the first dimension, that corresponds to the
% diagonal of a connectivity matrix (channel to itself interaction) in a
% vectorized compressed form.
% 
%       i = sc_bst_connect_vector_diagind(C)
%       i = sc_bst_connect_vector_diagind(C, isExpanded=0)
% 
% Input 'C' can either be a scalar (interpreted as the number of channels)
% or a vectorized compressed connectivity matrix of size NyxNtxNf, where Ny
% is the number of symetric channel combinations, from which the number of
% channel Nc is obtained from the formula Nc = (-1+sqrt(1+ 8*Ny))/2.
% 
% Optional parameter 'isExpanded' (default = 0) specifies if the vectorized
% form is expanded or not. If yes (1), Nc = sqrt(size(C,1)), if 'C' is not
% scalar.
% 
% Jonathan Godbout, 2013

if nargin<2, isExpanded = 0; end

if isscalar(x), Nc = x;
elseif ~isExpanded
    Nx = size(x,1);
    Nc = (-1+sqrt(1+ 8*Nx))/2;
elseif isExpanded
    Nc = sqrt(size(x,1));
end

if isExpanded
    i = find(eye(Nc)>0);
else
    a = 1:Nc;
    i = zeros(Nc,1);
    i(1) = 1;
    for k=2:Nc
        i(k) = i(k-1)+a(k);
    end
end

end
