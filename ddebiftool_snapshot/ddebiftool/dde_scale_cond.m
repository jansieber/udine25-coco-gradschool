function varargout=dde_scale_cond(varargin)
if isfield(varargin{1},'kind')
    varargout=cell(1,2);
    [varargout{:}]=loc_scale_cond(varargin{:});
else
    varargout{1}=dde_sys_cond_create('name',varargin{1},'reference',true,...
        'fun',@(p,pref)loc_scale_cond(p,pref,varargin{2:end}));
end
end
function [r,J]=loc_scale_cond(pt,pref,cond,varargin)
default={'res',1,'matrix',1,'ref',false};
options=dde_set_options(default,varargin,'pass_on');
[r0,J0]=condeval(cond,pt,pref);
if ~isempty(r0)
    if options.ref
        [rref,Jref]=condeval(cond,pref,pref);
        fac=1;
    else
        rref=r0;
        Jref=J0;
        fac=2;
    end
    prefac=rref'*options.matrix;
    r=prefac*r0-options.res;
    J=p_axpy(0,Jref(1),[]);
    for i=1:numel(J0)
        J=p_axpy(fac*prefac(i),Jref(i),J);
    end
else
    r=r0;
    J=J0;
end
end
%%
function [r,J]=condeval(cond,pt,pref)
if isstruct(cond)
    if cond.reference
        [r,J]=cond.fun(pt,pref);
    else
        [r,J]=cond.fun(pt);
    end
else
    [r,J]=cond(pt,pref);
end
end