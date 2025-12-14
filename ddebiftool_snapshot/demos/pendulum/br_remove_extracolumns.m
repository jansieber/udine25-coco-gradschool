function branch=br_remove_extracolumns(branch)
if ~isstruct(branch) || ~isfield(branch,'point')
    return
end
for i=1:length(branch.point)
    pt=branch.point(i);
    if isfield(pt,'nvec')&& isstruct(pt.nvec) && isfield(pt.nvec,'extracolumns')
        pt.nvec=rmfield(pt.nvec,'extracolumns');
    end
    branch.point(i)=pt;
end
end
