function pairname = sc_bst_connect_vector_label(name,rmdiag)
% Returns a vector of channel pairs labels of the form 'A-B'.  

if nargin<2, rmdiag = 0; end

Nc = numel(name);

channelpair = name(sc_bst_connect_index_ind2sub(Nc));
if rmdiag
	channelpair = sc_bst_connect_vector_rmdiag(channelpair);
end

Nx = size(channelpair,1);
pairname = cell(Nx,1);
for i=1:Nx
    pairname{i} = [channelpair{i,1},'-',channelpair{i,2}];
end

end