function c = sc_bst_process_subgroup(index)

if index > sc_bst_process_index(845)
    c = 'CARSM-B';
elseif index >= sc_bst_process_index(0)
    c = 'CARSM-A';
else
    c = 'CARSM-UNCLASSIFIED';
end

end