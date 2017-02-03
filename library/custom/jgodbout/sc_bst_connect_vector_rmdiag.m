function x = sc_bst_connect_vector_rmdiag(x)
% Removes diagonal elements (channel to itself interaction) from 
% connectivity matrix in a vectorized form.
% 
%       Co = sc_bst_connect_vector_rmdiag(Ci)
%       Co = sc_bst_connect_vector_rmdiag(Ci, isExpanded=0)
% 
% See <sc_bst_connect_vector_diagind.m> for input description.
% 
% Jonathan Godbout, 2013

x(sc_bst_connect_vector_diagind(x),:,:) = [];

end
