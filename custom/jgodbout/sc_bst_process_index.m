function index = sc_bst_process_index(index)
% Convert wanted process index (relative to the CARSM subcategory) to
% absolute process index (relative to all brainstorm index). Can be used to
% translate the process to the top or bottom of the list.

a = 1; % Scaling factor
b = -10000; % Translation factor
index = a*index+b;

end