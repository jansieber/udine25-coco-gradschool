function [kind,ptsol,userfuncs]=dde_get_kind(pt,funcs,arg)
%% get kind from point field or funcs
userfuncs=funcs;
if nargin<2 || isempty(funcs)|| ~isstruct(funcs) ||...
        ~isfield(funcs,'kind')
    kind=pt(1).kind;
    ptsol=pt;
else
    kind=funcs.kind;
    if nargout>1 && isfield(funcs,'get_comp')
        if nargin<3 
            ptsol=funcs.get_comp(pt,'solution');
        else
            ptsol=funcs.get_comp(pt,arg);
        end
    else
        ptsol=pt;
    end
    if nargout>2 && isfield(funcs,'userfuncs')
        userfuncs=funcs.userfuncs;
    end
end
end
