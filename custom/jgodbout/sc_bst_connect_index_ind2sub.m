function varargout = sc_bst_connect_index_ind2sub(Nc,ind,rmdiag)
% ONLY UPPER TRIANGLE

Nx = Nc*(Nc+1)/2;

if nargin<3, rmdiag = 0; end
if nargin<2, ind = 1:Nx; end
if islogical(ind), ind = find(ind); end

sqrind = find(triu(ones(Nc))>0);
Ni = min(Nx,numel(ind));
sub = zeros(Ni,2);
for k=1:Ni
    [sub(k,1),sub(k,2)] = ind2sub([Nc,Nc],sqrind(ind(k)));
end

if rmdiag
    if Ni~=Nx, error('Can not remove diagonal from incomplete vector'); end
    sub = sc_bst_connect_vector_rmdiag(sub);
end

if nargout<=1, varargout{1} = sub;
else
    varargout{1} = sub(:,1);
    varargout{2} = sub(:,2);
end

end
