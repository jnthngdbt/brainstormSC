BRAINSTORMSC, for "Brainstorm" add-on from "Sacré-Coeur" (Hôpital du Sacré-Coeur de Montréal, HSCM, more specifically the Center for Advanced Research in Sleep Medicine, CARSM, of HSCM). 

1. Make sure to have the latest version of Brainstorm installed (installed and in Matlab path). To install, download the toolbox from Brainstorm website. In the folder, there is a file (matlab function) named "brainstorm.m". Open it and run it in Matlab. It will install automatically (no need to add paths manually).

2. To install BRAINSTORMSC, run function "bst_sc_install.m" in its folder (do not add any path to Matlab, the function will do it automatically). It will copy the content of BRAINSTORMSC into Brainstorm folders (folder structure in BRAINTORMSC/library must match BRAINSTORM3). The functions that will be used are the ones copied in Brainstorm folders (not the ones in BRAINSTORMSC folders). This is an important fact only if you modify those functions.

3. For developpers (modifying BRAINSTORMSC functions in Brainstorm folders): run function "bst_sc_update.m" once in a while (if you modify the functions). It will update the functions of BRAINSTORMSC with their latest versions in Brainstorm folders (in the case they were modified). It is a 2-way update; newer versions overwrite older versions. Furthermore, new functions in BRAINSTORMSC folders will be added in Brainstorm folders.

Jonathan Godbout, M.Eng., 2013
