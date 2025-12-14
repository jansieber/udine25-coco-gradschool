function wrap=dde_psol_wrap(varargin)
%% wrap time points t into mesh periodically (including boundaries) 
% or map (possibly sparse) array y, which is supposed to have one column per entry in t,
% correspondingly. If y is given then t is assumed to be monotone
% increasing.
default={'mesh',zeros(0,1),'t',zeros(1,0),'itbd',NaN(1,2),'y',zeros(0,1)};
options=dde_set_options(default,varargin,'pass_on');
bdmsh=options.mesh([1,end]);
if ~isempty(options.t)
    wrap=min(max(mod(options.t-bdmsh(1),bdmsh(2)-bdmsh(1))+bdmsh(1),bdmsh(1)),bdmsh(2));
    wrap(options.t==bdmsh(2))=bdmsh(2);
    return
end
if ~isempty(options.y)
    itlen=options.itbd(2)-options.itbd(1)+1;
    [ir,ic,yvals]=find(options.y);
    ic=mod(ic-options.itbd(1)-1,itlen)+1;
    wrap=sparse(ir,ic,yvals,size(options.y,1),itlen);
    if ~issparse(options.y)
        wrap=full(wrap);
    end
end
end