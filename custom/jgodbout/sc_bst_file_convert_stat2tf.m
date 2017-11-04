function OutputFiles = sc_bst_file_convert_stat2tf( sInputs )
% sc_bst_file_convert_stat2tf: Converts "statistics" typed files to 
% "timefrequency" typed files. Useful until Brainstorm support connectivity
% statistics. The field TF contains statistics values (t-values or 
% p-values). Newly created files are add in database in same study has
% input files.
% 
% It has the format of a function instead of a process since (for now)
% process can't take statistics files.
%
% Authors: Jonathan Godbout, 2013

if ~strcmpi(unique({sInputs.FileType}),'ptimefreq')
    error('All files must be timefrequency statistics files');
end

OutputFiles = cell(numel(sInputs),1);
for iInput = 1:numel(sInputs)
    
    % Load statistics data
    sStat = in_bst_timefreq(sInputs(iInput).FileName);
    % Create time-frequency data from stat
    sTF = db_template('timefreqmat');
    for iField=intersect(fieldnames(sStat),fieldnames(sTF))'
        sTF.(iField{1}) = sStat.(iField{1});
    end
    sTF.TF          = sStat.tmap;
    sTF.DataType    = 'data';
    sTF.Options.statconvert = sStat;
%     sTF.Measure = 'other';
    sTF.Method = 'cohere';  % Otherwise display bug
    
    isConnectn = ~isempty(sTF.RefRowNames);
    if isConnectn,  FileType = 'timefreq_connectn';
    else            FileType = 'timefreq';
    end
    
    % ===== SAVE FILE =====
    % Study
    sStudy = bst_get('Study', sInputs(iInput).iStudy);
    % Output filename
    OutputFiles{iInput} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), FileType);
    % Save on disk
    bst_save(OutputFiles{iInput}, sTF, 'v6');
    % Register in database
    db_add_data(sInputs(iInput).iStudy, OutputFiles{iInput}, sTF);
end
    
% Update database interface
db_reload_studies([sInputs.iStudy]);

end
