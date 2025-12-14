function varargout=dde_extreme_cond(varargin)
%% provide condition for extrema of cx^(k)(t0)-d=0, with t0 a free parameter
% for periodic orbits (default k=0, change with option 'diff'
if isfield(varargin{1},'kind')
    varargout=cell(1,2);
    [varargout{:}]=loc_extreme_cond(varargin{:});
else
    varargout{1}=dde_sys_cond_create('name',varargin{1},'reference',false,...
        'fun',@(p)loc_extreme_cond(p,varargin{2:end}));
end
end
function [r,J]=loc_extreme_cond(p,c,ind_d, ind_t,varargin)
default={'diff',0};
options=dde_set_options(default,varargin,'pass_on');
t0=p.parameter(ind_t);
d=p.parameter(ind_d);
[x,Jx]=dde_coll_eva(p,mod(t0,1),'kron',true,'diff',options.diff);
[dx,dJx]=dde_coll_eva(p,mod(t0,1),'kron',true,'diff',options.diff+1);
d2x=dde_coll_eva(p,mod(t0,1),'kron',true,'diff',options.diff+2);
r=[c*x-d;c*dx];
J=repmat(p_axpy(0,p,[]),2,1);
J(1).profile=reshape(c*Jx,size(p.profile));
J(1).parameter(ind_d)=-1;
J(1).parameter(ind_t)=c*dx;
J(2).profile=reshape(c*dJx,size(p.profile));
J(2).parameter(ind_t)=c*d2x;
end

