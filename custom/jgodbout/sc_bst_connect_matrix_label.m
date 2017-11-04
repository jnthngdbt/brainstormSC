function pairname = sc_bst_connect_matrix_label(name)
% Returns a vector of channel pairs labels of the form 'A-B'. 

Nc = numel(name);

pairname = cell(Nc,Nc);
for i=1:Nc
    for j=1:Nc
        pairname{i,j} = sprintf('%s-%s',name{i},name{j});
    end
end

end