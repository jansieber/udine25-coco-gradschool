%% residual & derivative wrt locating extrema of functions along coll profile
% function has form F(tz)=f(x(tz),x(tz-tau_1)...,p),
% and condition is F(tz)=F'(tz)=0 with additional unknown tz,
% formulated in sys_cond format cond(x(tz),x(tz-tau1),...,x(tz-taud),p)=0
% for some t0 and [d/dt] cond(x(t),x(t-tau1),...,x(t-taud),p)|_{t=t0}=0
% This requires the additional parameter tz ,which is assumed to be already
% included in point.parameter at position 'itz' (optional argument, default
% is lasst parameter)
%
% inputs
%
% * funcs: system functions
% * point: point (type psol) for cond(x,p)=0 is to be determined
% * cond: scalar function c(xx,p) where xx is n x (d+1), p are
% parameters (d=number of delays)
% 
% Optional/named  inputs
%
% * additional args for set_funcs
% *'itz'=position of parameter index for additional parameter t0 (default=end of parameter
% list for p)
% * 'dcond': directional derivative of cond:dcond(x,p,dx,dp), if empty,
% finite difference is used
% * 'freepar': Jpt only contains derivative wrt freepar, default if not
% provided: all parameters
%
% outputs
% r residual: F(tz) and F'(tz) for psol
% J Jacobian (struct of type point): (2x1) point structure for psol
%%
function varargout=dde_nlin_extreme_cond(varargin)
if isfield(varargin{1},'kind')
    varargout=cell(1,2);
    [varargout{:}]=loc_nlin_extreme_cond(varargin{:});
else
    varargout{1}=dde_sys_cond_create('name',varargin{1},'reference',false,...
        'fun',@(p)loc_nlin_extreme_cond(p,varargin{2:end}));
end
end
function [r,Jpt]=loc_nlin_extreme_cond(point,funcs,cond,varargin)
ntau=dde_num_delays(funcs);
npar=length(point.parameter);
default={'itz',npar,'free_par',1:npar,'ntau',ntau,'dcond',[],'diff',0};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
nx=size(point.profile,1);
ptext=point;
ptext.profile=[point.profile;zeros(1,length(point.mesh))];
ptext.parameter(options.itz)=[];
free_par_ext=options.free_par;
free_par_ext(options.free_par==options.itz)=[];
sys_tau=@(it,x,p)funcs.wrap_tau(it,x(1:end-1,:,:),p);
sys_dirdtau=arrayfun(@(k){@(it,xx,p,dxx,dp)funcs.wrap_dirdtau{k}(it,...
    xx(1:end-1,:,:),p,dxx(1:end-1,:,:),dp)},1:length(funcs.wrap_dirdtau));
sys_ntau=funcs.sys_ntau;
if isempty(options.dcond)
    deriarg={};
else
    deriarg={'sys_dirderi',{@(x,p,dx,dp)loc_dirderi(x,p,dx,dp,options.dcond)}};
end
tfuncs=set_funcs('sys_rhs',@(x,p)loc_rhs(x,p,cond),deriarg{:},...
    'wrap_tau',sys_tau,'sys_ntau',sys_ntau,'wrap_dirdtau',sys_dirdtau,...
    'sys_tau_seq',funcs.sys_tau_seq,...
    'lhs_matrix',zeros(1,nx+1),...
    'x_vectorized',funcs.x_vectorized,'p_vectorized',funcs.p_vectorized,pass_on{:});
mth=getfield(df_mthod('psol'),'point');
mth.matrix='sparse';
[Jt,rt]=dde_coll_jac(tfuncs,ptext,free_par_ext,...
    'collocation_parameters',ptext.mesh,'c_is_tvals',true);
cprofile=rt.';
collargs={setfield(point,'profile',cprofile),point.parameter(options.itz)}; %#ok<*SFLD>
r=NaN(3,1);
for i=1:3
    r(i)=dde_coll_eva(collargs{:},'diff',i-1+options.diff);
end
r2=r(3);
r=r(1:2);
if nargout==1
    return
end
ind=dde_ind_from_point(ptext,free_par_ext);
Jx=Jt.profile(:,ind.profile(1:nx,:));
Jdef=Jt.profile(:,ind.profile(nx+1,:));
Jpt=repmat(p_axpy(0,point,[]),2,1);
jtzargs={setfield(point,'profile',cprofile),point.parameter(options.itz),'output','matrix'};
for i=1:2
    Jtz=dde_coll_eva(jtzargs{:},'diff',i-1+options.diff);
    Jtz_inv=-Jtz/Jdef;
    Jpt(i).profile(:)=Jtz_inv*Jx;
    Jpt(i).parameter(free_par_ext)=Jtz_inv*Jt.parameter;
    Jpt(i).period=Jtz_inv*Jt.period;
end
Jpt(1).parameter(options.itz)=r(2);
Jpt(2).parameter(options.itz)=r2;
end
%%
function r=loc_rhs(xx,p,cond)
nvec=size(xx,3);
def=reshape(xx(end,1,:),1,nvec);
x=xx(1:end-1,:,:);
r=cond(x,p)-def;
end
%%
function r=loc_dirderi(xx,p,dxx,dp,dcond)
% sysdtau=@(itau,x,p,nx,np)
nvec=size(xx,3);
ddef=reshape(dxx(end,1,:),1,nvec);
x=xx(1:end-1,:,:);
dx=dxx(1:end-1,:,:);
r=dcond(x,p,dx,dp)-ddef;
end
