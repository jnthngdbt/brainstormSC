function sc_bst_update
% Function that synchronizes 'Brainstorm SC Toolbox' functions that are in
% 'Brainstorm' folders (the ones that are actually used, in Matlab path)
% with their equivalent in 'Brainstorm SC Toolbox' folders. 
% 
% It is a 2-way update, so newer files overwrites older files.
% 
% It only calls function SC_BST_INSTALL (run when installed the toolbox),
% which does exactly what we want here. A distinct function 'update' avoids
% confusion and invite user to run this function as often as possible when
% developping.
% 
% IMPORTANT NOTE: RUN AS OFTEN AS POSSIBLE, SINCE WHEN UPDATING BRAINSTORM,
% IT DELETES ANY OTHER FILES THAT WERE ADDED OR MODIFIED. THE 'Brainstorm
% SC Toolbox' MUST THEN BE CONSTANTLY UPDATED TO KEEP LATEST VERSIONS OF
% WORKING FUNCTIONS (LOCATED IN BRAINSTORM PATH).
% 
% Jonathan Godbout, 2013

sc_bst_install;

end