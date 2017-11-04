function sc_bst_connect_stat_table(s,pmax)
% Prints statistics values in command window. Input S is a struct manually
% exported from Brainstorm (right-click on stat data node, then File >
% Export to Matlab)

if ~isfield(s,'RowNames') || ~isfield(s,'pmap')
    error('Input must be connectivity statistics');
end
if numel(s.RowNames)~=numel(s.RefRowNames)
    error('Unsupported 1xN connectivity yet');
end

if nargin<2, pmax = 0.9; end
% if nargin<2, sortby = 'pmap'; end

[Nx,Nt,Nf] = size(s.pmap);

labelPairs = sc_bst_connect_vector_label(s.RowNames);
labelPairs = repmat(labelPairs(:),Nf,1);

freqBand = cell(1,1,Nf);
freqBand(:) = s.Freqs(:,2);
freqBand = repmat(freqBand,[Nx,1,1]);
freqBand = freqBand(:);

pmap = s.pmap(:);
tmap = s.tmap(:);
df = s.df(:);
mask = s.Options.statcond.mask(:);
mask_mp = s.Options.statcond.mask_mp(:);
alpha = s.Options.statcond.alpha;
mp = s.Options.statcond.mpmethod;

tablehdr = {'channels','p-values',    'score',       sprintf('p<%.3f',alpha),sprintf('p<%.3f (%s)',alpha,mp),'DoF',       'Band (Hz)'};
table    = [labelPairs,num2cell(pmap),num2cell(tmap),num2cell(mask),         num2cell(mask_mp),              num2cell(df),freqBand];

table = table(pmap<=pmax,:);

[sortval,sortidx] = sort(cell2mat(table(:,2)),'ascend');
table = table(sortidx,:);

disp(s.Comment)
disp([tablehdr;table])

% % Display
% tablestr = [];
% tablestr = [tablestr; sprintf('%15s,')];
% 
% tablestr

end
