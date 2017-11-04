function m = sc_bst_connect_vector_mask_pair(label,pair,isSym)
% Return a logical mask M of channel pairs PAIR. PAIR is a cell array of
% wanted channel pairs labels (ex: 'F3-F4'). ISSYM is a flag of pair
% symmetry (if =1, 'F3-F4'=='F4-F3'). LABEL is a cell array of channels
% labels, needed to know the channels and pairs order.

if nargin<3, isSym = 1; end

if ischar(pair), pair = {pair}; end

if isempty(pair), pair = sc_bst_connect_vector_label(label); end

PairMatrix = sc_bst_connect_matrix_label(label);
iPairFilter = false(size(PairMatrix));
for i=1:numel(pair)
    iPairFilter = iPairFilter | strcmpi(PairMatrix, pair{i});
end

if isSym
    iPairFilter = iPairFilter | iPairFilter'; % Undirected connectivity
end

m = sc_bst_connect_format_mat2vec(iPairFilter)>0;

end