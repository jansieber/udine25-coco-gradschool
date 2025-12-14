function varargout=POEV1Init(funcs,point,method,varargin)
%% crude initial guess for fold of periodic orbits
%
%
ip=funcs.ip;
default={'v_scal',@(p,pref)sys_cond_POEV1_norm(p,ip,[ip.periodnullpar,ip.nullparind(:,2)'],'res',0),'nulldim',1};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
out=cell(1,3);
[out{1:3}]=svd_coll_init(funcs,ip,point,0,options.nulldim,method.point,'v_scal',options.v_scal,pass_on{:});
varargout=out;
end