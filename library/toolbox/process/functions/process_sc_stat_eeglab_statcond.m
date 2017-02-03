function varargout = process_sc_stat_eeglab_statcond( varargin )
% process_sc_stat_eeglab_statcond: Calls EEGLAB function STATCOND
%
% Authors: Jonathan Godbout, 2013

macro_methodcall;
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'STATISTICS > EEGLAB > Conditions comparisons statistics (statcond.m)';
    sProcess.FileTag     = '';
    sProcess.Category    = 'Custom';
    sProcess.Index       = sc_bst_process_index(1545);
    sProcess.SubGroup    = sc_bst_process_subgroup(sProcess.Index);
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    %
    sProcess.options.method.Comment = 'Method for computing p-values: ';
    sProcess.options.method.Type    = 'combobox';
    sProcess.options.method.Value   = {3,{'parametric','permutation','bootstrap'}};
    %
    sProcess.options.tail.Comment = '1 tail or 2 tails: ';
    sProcess.options.tail.Type    = 'combobox';
    sProcess.options.tail.Value   = {2,{'one','both','upper','lower'}};
    %
    sProcess.options.mpmethod.Comment = 'Method for multiple comparisons correction: ';
    sProcess.options.mpmethod.Type    = 'combobox';
    sProcess.options.mpmethod.Value   = {1,{'FDR','Bonferroni','Holm-Bonferroni','Hochberg','None'}};
    %
    sProcess.options.variance.Comment = 'Variance estimation (for degree of freedom): ';
    sProcess.options.variance.Type    = 'combobox';
    sProcess.options.variance.Value   = {1,{'inhomogenous','homogenous'}};
    %
    sProcess.options.forceunpaired.Comment = 'Force un-pair data (ignored for non-equal sized data)';
    sProcess.options.forceunpaired.Type    = 'checkbox';
    sProcess.options.forceunpaired.Value   = 0;
    %
    sProcess.options.forceanova.Comment = 'Force the use of ANOVA calculation even for 2x1 designs';
    sProcess.options.forceanova.Type    = 'checkbox';
    sProcess.options.forceanova.Value   = 0;
    %
    sProcess.options.alpha.Comment = 'p-value threshold value';
    sProcess.options.alpha.Type    = 'value';
    sProcess.options.alpha.Value   = {0.05,'',5};
    %
    sProcess.options.naccu.Comment = ' Number of surrogate data copies to use in "perm" or "bootstrap" method estimation';
    sProcess.options.naccu.Type    = 'value';
    sProcess.options.naccu.Value   = {500,'',0};
    % Specific frequency bands
    sProcess.options.freqfilter.Comment = 'Only consider specific frequency bands (All=0)';
    sProcess.options.freqfilter.Type    = 'value';
    sProcess.options.freqfilter.Value   = {[1 2],'',0};
    % Specific sensors
    sProcess.options.chanfilter.Comment = 'Only consider specific sensors (ex: ''F3,F4'' or ''F3-P3,F3-F4'' for connectivity).';
    sProcess.options.chanfilter.Type    = 'text';
    sProcess.options.chanfilter.Value   = [];
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>   

OutputFiles = [];

% Conditions
cond = {sInputs.Condition};
condNames = unique(cond);
nCond = numel(condNames);

if nCond==1
    strMsg = 'Only 1 condition, can not compute comparison stats!';
    bst_report('Error',sProcess, sInputs, strMsg);
    return;
end

% Specific frequency band index
iFreq = sProcess.options.freqfilter.Value{1};
iFreq(iFreq==0) = [];

% Specific channels
ChanFilter = sProcess.options.chanfilter.Value;
if ~isempty(ChanFilter)
    ChanFilter = textscan(ChanFilter,'%s','delimiter',',');
    ChanFilter = ChanFilter{1};
end

% Extract condition specific data (stat along 4th dimension)
data = cell(1,nCond);
for iInput=1:numel(sInputs)
    [sMat,nameMat] = in_bst(sInputs(iInput).FileName);
    
    iCond = strcmpi(condNames, sInputs(iInput).Condition);    
    datai = sMat.(nameMat); 
        
   data{iCond} = cat(4,data{iCond}, datai);
    
end
[Nc,Nt,Nf] = size(datai);

% Unique file names (reflects preprocessing done)
fileComment = unique({sInputs.Comment});
% if numel(fileComment)>1
%     bst_report('Error', sProcess, sInputs, 'Non-unique file comments; may reflect different pre-processing hence different things to statisticaly compare.')
% end

% Filetype
isData = 0;
isResults = 0;
isTimefreq = 0;
isConnectivity = 0;
% Initialize output structure
sOutput = db_template('statmat');
% ALL INPUT FILES MUST BE OF SAME TYPE WITH SAME CHANNEL-TIME-FREQ DIMS
switch lower(sInputs(1).FileType)
    case 'data'
        isData = 1;
        DataMat = in_bst_data(sInputs(1).FileName, 0, 'Time');
        sOutput.ColormapType = 'stat2';
        sOutput.Time = DataMat.Time;
    case 'results'
        isResults = 1;
        sOutput.ColormapType = 'stat1';
        % Read some info from the first file
        ResultsMat = in_bst_results(sInputs(1).FileName, 0, 'Atlas', 'SurfaceFile','Time');
        sOutput.Atlas       = ResultsMat.Atlas;
        sOutput.SurfaceFile = ResultsMat.SurfaceFile;
        sOutput.Time        = ResultsMat.Time;
    case 'timefreq'
        sOutput.ColormapType = 'stat1';
        % Read some info from the first file
        TimefreqMat = in_bst_timefreq(sInputs(1).FileName, 0, 'Atlas', 'DataType', 'SurfaceFile', 'TimeBands', 'Freqs', 'RefRowNames', 'RowNames', 'Measure', 'Method', 'Options','Time','DataFile');
        sOutput.Atlas       = TimefreqMat.Atlas;
        sOutput.DataType    = TimefreqMat.DataType;
        sOutput.SurfaceFile = TimefreqMat.SurfaceFile;
        sOutput.TimeBands   = TimefreqMat.TimeBands;
        sOutput.Freqs       = TimefreqMat.Freqs;
        sOutput.RefRowNames = TimefreqMat.RefRowNames;
        sOutput.RowNames    = TimefreqMat.RowNames;
        sOutput.Measure     = 'other';
        sOutput.Method      = 'ttest';
        sOutput.Options     = TimefreqMat.Options;
        sOutput.Time        = TimefreqMat.Time;   
        
        if isempty(sOutput.RefRowNames)
            isTimefreq = 1;
        else
            isConnectivity = 1;
        end        
end

% Build frequency band mask
maskFreq = true(Nc,Nt,Nf);
if (isTimefreq || isConnectivity) && ~isempty(iFreq)
    maskFreq = false(Nc,Nt,Nf);
    maskFreq(:,:,iFreq) = true;
end

% Build channel (or pairs for connectivity) mask
maskChannel = true(Nc,Nt,Nf);
if ~isempty(ChanFilter)
    if isConnectivity
        maskChannel = sc_bst_connect_vector_mask_pair(sOutput.RowNames,ChanFilter,1);
        maskChannel = repmat(maskChannel(:),[1,Nt,Nf]);
    elseif isTimefreq
        maskChannel = cellfun(@(x)ismember(x,upper(ChanFilter)),upper(sOutput.RowNames));
        maskChannel = repmat(maskChannel(:),[1,Nt,Nf]);
    end
end

% % === COMPUTE TEST ===

% Get parameters
mpmethod = sProcess.options.mpmethod.Value{2}{sProcess.options.mpmethod.Value{1}}
method = sProcess.options.method.Value{2}{sProcess.options.method.Value{1}}
tail = sProcess.options.tail.Value{2}{sProcess.options.tail.Value{1}}
variance = sProcess.options.variance.Value{2}{sProcess.options.variance.Value{1}}
naccu  = sProcess.options.naccu.Value{1}
alpha  = sProcess.options.alpha.Value{1}
paired = sProcess.options.forceunpaired.Value==0
if paired, paired = 'on'; else paired = 'off'; end
forceanova = sProcess.options.forceanova.Value
if forceanova, forceanova = 'on'; else forceanova = 'off'; end
% tail = 'one';

% Compute using EEGLAB function. All tests possible.
% [stats, df, pvals, surrog] = statcond( data,'mode',method,'paired',paired,'naccu',naccu,'tail','one');
S = statcond( data,'mode',method,'paired',paired,'naccu',naccu,'tail',tail,'forceanova',forceanova,'alpha',alpha,'variance',variance,'structoutput','on');

% Compute filter mask. All p-values not in mask will be set to 0, but the
% stat will keep their values
maskNull = isnan(S.stat); % Null variances
FilterMask = maskChannel & maskFreq & ~maskNull;
S.stat(maskNull) = 0;
S.pval(maskNull) = 1;
S.pval(~FilterMask) = 1;

% Determine which analysis was done
isPaired = strcmpi(paired,'on') && numel(unique(cellfun(@(x)size(x,3), data)))==1;
if isPaired && isfield(S,'t') && size(data,1)==1 && size(data,2)==2
    comment = sprintf('Paired t-test [%s]-[%s]',condNames{:});
elseif ~isPaired && isfield(S,'t') && size(data,1)==1 && size(data,2)==2
    comment = sprintf('Unpaired t-test [%s]-[%s]',condNames{:});
elseif isPaired && isfield(S,'f') && size(data,1)==1 && size(data,2)>2
    comment = 'Paired 1-way ANOVA';
elseif ~isPaired && isfield(S,'f') && size(data,1)==1 && size(data,2)>2
    comment = 'Unpaired 1-way ANOVA';
else
    bst_report('Error',sProcess, sInputs, 'Unsupported test');
    return;
end
sProcess.Comment = [comment];
% sProcess.Comment = [comment,' | ',fileComment{1}];

% Correct for multiple comparisons (Nan p-values are not considered)
mask_mp = false(size(S.pval));
mask_mp(FilterMask) = MultipleComparisons(S.pval(FilterMask),alpha,mpmethod);

% === FILL RESULT STRUCTURE ===
sOutput.df   	= S.df;
sOutput.tmap   	= S.stat;
sOutput.pmap   	= S.pval;
sOutput.Options.statcond.mpmethod = mpmethod;
sOutput.Options.statcond.method = method;
sOutput.Options.statcond.naccu = naccu;
sOutput.Options.statcond.alpha = alpha;
sOutput.Options.statcond.paired   	= paired;
sOutput.Options.statcond.tail   	= tail;
sOutput.Options.statcond.mask   	= S.mask;
sOutput.Options.statcond.mask_mp = mask_mp;
sOutput.Options.statcond.nDiscoveries = numel(find(mask_mp));
sOutput.Options.statcond.ci   	= S.ci;
sOutput.Options.statcond.filtermask   	= FilterMask;
sOutput.Options.statcond.channelsubmask   	= maskChannel;
sOutput.Options.statcond.freqsubmask   	= maskFreq;

% SAVE (code extracted from BST_PROCESS>STAT1, since connectivity
% statistics not supported yet (8-july-2013)
OutputFiles = SaveProcessStat(sProcess, sInputs, [], sOutput);

% % % % % % Connectivity statistics are not supported yet, so create a dummy
% % % % % % connectivity data containing p-values instead of connectivity values (in
% % % % % % order to visualize p-values using built-in connectivity visualization
% % % % % if exist('TimefreqMat') % Assume it is connectivity
% % % % %     SaveConnectStatForVis(sProcess, sInputs, [], sOutput,TimefreqMat);
% % % % % end

end

%% ===== PROCESS: STAT =====
function OutputFiles = SaveProcessStat(sProcess, sInputA, sInputB, sOutput)

% % % % % % % %     % ===== CALL PROCESS =====
% % % % % % % %     isStat1 = strcmpi(sProcess.Category, 'Stat1');
% % % % % % % %     if isStat1
% % % % % % % %         sOutput = sProcess.Function('NoCatch', 'Run', sProcess, sInputA);
% % % % % % % %     else
% % % % % % % %         sOutput = sProcess.Function('NoCatch', 'Run', sProcess, sInputA, sInputB);
% % % % % % % %     end
% % % % % % % %     if isempty(sOutput)
% % % % % % % %         OutputFiles = {};
% % % % % % % %         return;
% % % % % % % %     end
    
    isStat1 = 1;

    % ===== GET OUTPUT STUDY =====
    % Display waitbar
    bst_progress('text', 'Saving results...');
    % Get number of subjects that are involved
    if isStat1
        uniqueSubjectName = unique({sInputA.SubjectFile});
        uniqueStudy       = unique([sInputA.iStudy]);
    else
        uniqueSubjectName = unique([{sInputA.SubjectFile}, {sInputB.SubjectFile}]);
        uniqueStudy       = unique([sInputA.iStudy, sInputB.iStudy]);
    end
    % If all files share same study: save in it
    if (length(uniqueStudy) == 1)
        [sStudy, iStudy] = bst_get('Study', uniqueStudy);
    % If all files share the same subject: save in intra-analysis
    elseif (length(uniqueSubjectName) == 1)
        % Get subject
        [sSubject, iSubject] = bst_get('Subject', uniqueSubjectName{1});
        % Get intra-subjet analysis study for this subject
        [sStudy, iStudy] = bst_get('AnalysisIntraStudy', iSubject);
    else
        % Get inter-subjet analysis for this subject
        [sStudy, iStudy] = bst_get('AnalysisInterStudy');
    end

    % ===== CREATE OUTPUT STRUCTURE =====
    % Template structure for stat files
    sOutput.Type = sInputA(1).FileType;
    % Comment: forced in the options
    if isfield(sProcess.options, 'Comment') && isfield(sProcess.options.Comment, 'Value') && ~isempty(sProcess.options.Comment.Value)
        sOutput.Comment = sProcess.options.Comment.Value;
    % Regular comment
    else
        sOutput.Comment = sProcess.Function('FormatComment', sProcess);
        if ~isStat1
            % Get comment for files A and B
            [tmp__, tmp__, CommentA] = bst_process('GetOutputStudy', sProcess, sInputA);
            [tmp__, tmp__, CommentB] = bst_process('GetOutputStudy', sProcess, sInputB);
            % Get full comment
            sOutput.Comment = [sOutput.Comment ': ' CommentA ' vs. ' CommentB];
        end
    end
    % Results: Get extra infotmation
    if strcmpi(sInputA(1).FileType, 'results')
        % Load extra fields
        ResultsMat = in_bst_results(sInputA(1).FileName, 0, 'HeadModelType', 'SurfaceFile', 'nComponents');
        % Copy fields
        sOutput.HeadModelType = ResultsMat.HeadModelType;
        sOutput.SurfaceFile   = ResultsMat.SurfaceFile;
        sOutput.nComponents   = ResultsMat.nComponents;
    end
    % History
    sOutput = bst_history('add', sOutput, 'stat', sProcess.Comment);
    % History: List files A
    sOutput = bst_history('add', sOutput, 'stat', 'List of files in group A:');
    for i = 1:length(sInputA)
        sOutput = bst_history('add', sOutput, 'stat', [' - ' sInputA(i).FileName]);
    end
    % History: List files B
    sOutput = bst_history('add', sOutput, 'stat', 'List of files in group B:');
    for i = 1:length(sInputB)
        sOutput = bst_history('add', sOutput, 'stat', [' - ' sInputB(i).FileName]);
    end
    
    % ===== SAVE FILE =====
    % Output filetype
    fileType = ['p', bst_process('GetFileTag', sInputA(1).FileName)];
    % Output filename
    OutputFiles{1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), fileType);
    % Save on disk
    bst_save(OutputFiles{1}, sOutput, 'v6');
    % Register in database
    db_add_data(iStudy, OutputFiles{1}, sOutput);
    
end

function [mask] = MultipleComparisons(pval,alpha,method)
% Adapted from Fieldtrip's FT_STATISTICS_MONTECARLO

stat.prob = pval;
cfg.alpha = alpha;

switch lower(method)
%   case 'max'
%     % the correction is implicit in the method
%     fprintf('using a maximum-statistic based method for multiple comparison correction\n');
%     fprintf('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
%     stat.mask = stat.prob<=cfg.alpha;
%     stat.posdistribution = posdistribution;
%     stat.negdistribution = negdistribution;
%   case 'cluster'
%     % the correction is implicit in the method
%     fprintf('using a cluster-based method for multiple comparison correction\n');
%     fprintf('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
%     stat.mask = stat.prob<=cfg.alpha;
  case 'bonferroni'
    fprintf('performing Bonferroni correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    stat.mask = stat.prob<=(cfg.alpha ./ numel(stat.prob));
  case {'holm','holm-bonferroni'}
    % test the most significatt significance probability against alpha/N, the second largest against alpha/(N-1), etc.
    fprintf('performing Holm-Bonferroni correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    [pvals, indx] = sort(stat.prob(:));                                   % this sorts the significance probabilities from smallest to largest
    k = find(pvals > (cfg.alpha ./ ((length(pvals):-1:1)')), 1, 'first'); % compare each significance probability against its individual threshold
    mask = (1:length(pvals))'<k;   
    stat.mask = zeros(size(stat.prob));
    stat.mask(indx) = mask;
  case 'hochberg'
    % test the most significatt significance probability against alpha/N, the second largest against alpha/(N-1), etc.
    fprintf('performing Hochberg''s correction for multiple comparisons (this is *not* the Benjamini-Hochberg FDR procedure!)\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    [pvals, indx] = sort(stat.prob(:));                     % this sorts the significance probabilities from smallest to largest
    k = find(pvals <= (cfg.alpha ./ ((length(pvals):-1:1)')), 1, 'last'); % compare each significance probability against its individual threshold
    mask = (1:length(pvals))'<=k;   
    stat.mask = zeros(size(stat.prob));
    stat.mask(indx) = mask;    
  case 'fdr'
    fprintf('performing FDR correction for multiple comparisons\n');
    fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
    stat.mask = fdr(stat.prob, cfg.alpha);
  otherwise
    fprintf('not performing a correction for multiple comparisons\n');
    stat.mask = stat.prob<=cfg.alpha;
end

mask = stat.mask;

end


function [h] = fdr(p, q)

% FDR false discovery rate
%
% Use as
%   h = fdr(p, q)
%
% This implements
%   Genovese CR, Lazar NA, Nichols T.
%   Thresholding of statistical maps in functional neuroimaging using the false discovery rate.
%   Neuroimage. 2002 Apr;15(4):870-8.

% Copyright (C) 2005, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: fdr.m 7123 2012-12-06 21:21:38Z roboos $

% convert the input into a row vector
dim = size(p);
p = reshape(p, 1, numel(p));

% sort the observed uncorrected probabilities
[ps, indx] = sort(p);

% count the number of voxels
V = length(p);

% compute the threshold probability for each voxel
pi = ((1:V)/V)  * q / c(V);

h = (ps<=pi);

% undo the sorting
[dum, unsort] = sort(indx);
h = h(unsort);

% convert the output back into the original format
h = reshape(h, dim);

end

function s = c(V)
% See Genovese, Lazar and Holmes (2002) page 872, second column, first paragraph
if V<1000
  % compute it exactly
  s = sum(1./(1:V));
else
  % approximate it
  s = log(V) + 0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951463144724980708248096050401448654283622417399764492353625350033374293733773767394279259525824709491600873520394816567;
end

end

% % % % % % ========================================================================
% % % % % function SaveConnectStatForVis(sProcess, sInputA, sInputB, sOutput,TimefreqMat)
% % % % % 
% % % % %     TimefreqMat.TF = sOutput.pmap;
% % % % %     sOutput = TimefreqMat;
% % % % %     sOutput.Comment = 'DUMMY CONNECTIVITY FOR P-VALUES VISUALIZATION';
% % % % %     sOutput.ColormapType = 'stat1';
% % % % % 
% % % % %     % ===== GET OUTPUT STUDY =====
% % % % %     % Display waitbar
% % % % %     bst_progress('text', 'Saving results...');
% % % % %     % Get number of subjects that are involved
% % % % %         uniqueSubjectName = unique({sInputA.SubjectFile});
% % % % %         uniqueStudy       = unique([sInputA.iStudy]);
% % % % %     % If all files share same study: save in it
% % % % %     if (length(uniqueStudy) == 1)
% % % % %         [sStudy, iStudy] = bst_get('Study', uniqueStudy);
% % % % %     % If all files share the same subject: save in intra-analysis
% % % % %     elseif (length(uniqueSubjectName) == 1)
% % % % %         % Get subject
% % % % %         [sSubject, iSubject] = bst_get('Subject', uniqueSubjectName{1});
% % % % %         % Get intra-subjet analysis study for this subject
% % % % %         [sStudy, iStudy] = bst_get('AnalysisIntraStudy', iSubject);
% % % % %     else
% % % % %         % Get inter-subjet analysis for this subject
% % % % %         [sStudy, iStudy] = bst_get('AnalysisInterStudy');
% % % % %     end
% % % % % 
% % % % %     % ===== CREATE OUTPUT STRUCTURE =====
% % % % % % %     % Template structure for stat files
% % % % % % %     sOutput.Type = sInputA(1).FileType;
% % % % % % %     % Comment: forced in the options
% % % % % % %     if isfield(sProcess.options, 'Comment') && isfield(sProcess.options.Comment, 'Value') && ~isempty(sProcess.options.Comment.Value)
% % % % % % %         sOutput.Comment = sProcess.options.Comment.Value;
% % % % % % %     % Regular comment
% % % % % % %     else
% % % % % % %         sOutput.Comment = sProcess.Function('FormatComment', sProcess);
% % % % % % %     end
% % % % % % %     % Results: Get extra infotmation
% % % % % % %     if strcmpi(sInputA(1).FileType, 'results')
% % % % % % %         % Load extra fields
% % % % % % %         ResultsMat = in_bst_results(sInputA(1).FileName, 0, 'HeadModelType', 'SurfaceFile', 'nComponents');
% % % % % % %         % Copy fields
% % % % % % %         sOutput.HeadModelType = ResultsMat.HeadModelType;
% % % % % % %         sOutput.SurfaceFile   = ResultsMat.SurfaceFile;
% % % % % % %         sOutput.nComponents   = ResultsMat.nComponents;
% % % % % % %     end
% % % % % 
% % % % %     % ===== SAVE FILE =====
% % % % %     % Output filetype
% % % % % % % % %     fileType = ['p', bst_process('GetFileTag', sInputA(1).FileName)];
% % % % %     fileType = [bst_process('GetFileTag', sInputA(1).FileName)];
% % % % %     % Output filename
% % % % %     OutputFiles{1} = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), fileType);
% % % % %     % Save on disk
% % % % %     bst_save(OutputFiles{1}, sOutput, 'v6');
% % % % %     % Register in database
% % % % %     db_add_data(iStudy, OutputFiles{1}, sOutput);
% % % % %     
% % % % % end


% % % % % FOLLOWING CODE EXTRACTED FROM FIELDTRIP'S FT_STATISTICS_MONTECARLO
% % % % switch lower(cfg.correctm)
% % % %   case 'max'
% % % %     % the correction is implicit in the method
% % % %     fprintf('using a maximum-statistic based method for multiple comparison correction\n');
% % % %     fprintf('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
% % % %     stat.mask = stat.prob<=cfg.alpha;
% % % %     stat.posdistribution = posdistribution;
% % % %     stat.negdistribution = negdistribution;
% % % %   case 'cluster'
% % % %     % the correction is implicit in the method
% % % %     fprintf('using a cluster-based method for multiple comparison correction\n');
% % % %     fprintf('the returned probabilities and the thresholded mask are corrected for multiple comparisons\n');
% % % %     stat.mask = stat.prob<=cfg.alpha;
% % % %   case 'bonferroni'
% % % %     fprintf('performing Bonferroni correction for multiple comparisons\n');
% % % %     fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
% % % %     stat.mask = stat.prob<=(cfg.alpha ./ numel(stat.prob));
% % % %   case 'holm'
% % % %     % test the most significatt significance probability against alpha/N, the second largest against alpha/(N-1), etc.
% % % %     fprintf('performing Holm-Bonferroni correction for multiple comparisons\n');
% % % %     fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
% % % %     [pvals, indx] = sort(stat.prob(:));                                   % this sorts the significance probabilities from smallest to largest
% % % %     k = find(pvals > (cfg.alpha ./ ((length(pvals):-1:1)')), 1, 'first'); % compare each significance probability against its individual threshold
% % % %     mask = (1:length(pvals))'<k;   
% % % %     stat.mask = zeros(size(stat.prob));
% % % %     stat.mask(indx) = mask;
% % % %   case 'hochberg'
% % % %     % test the most significatt significance probability against alpha/N, the second largest against alpha/(N-1), etc.
% % % %     fprintf('performing Hochberg''s correction for multiple comparisons (this is *not* the Benjamini-Hochberg FDR procedure!)\n');
% % % %     fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
% % % %     [pvals, indx] = sort(stat.prob(:));                     % this sorts the significance probabilities from smallest to largest
% % % %     k = find(pvals <= (cfg.alpha ./ ((length(pvals):-1:1)')), 1, 'last'); % compare each significance probability against its individual threshold
% % % %     mask = (1:length(pvals))'<=k;   
% % % %     stat.mask = zeros(size(stat.prob));
% % % %     stat.mask(indx) = mask;    
% % % %   case 'fdr'
% % % %     fprintf('performing FDR correction for multiple comparisons\n');
% % % %     fprintf('the returned probabilities are uncorrected, the thresholded mask is corrected\n');
% % % %     stat.mask = fdr(stat.prob, cfg.alpha);
% % % %   otherwise
% % % %     fprintf('not performing a correction for multiple comparisons\n');
% % % %     stat.mask = stat.prob<=cfg.alpha;
% % % % end




