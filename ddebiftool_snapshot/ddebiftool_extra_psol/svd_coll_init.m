function varargout=svd_coll_init(funcs,ip,point,omega,nulldim,method,varargin)
%% nullspace of linearized system in psol
%
%
internalpar=ip.nullparind(:,1)';
nint=length(internalpar);
default={'v_scal',@(p,pref)sys_cond_lincoll_norm(p,ip,'res',0),'output','ini'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% iniitalize extended point
ptini=dde_collbif_template(ip,point,omega);
%% which parameters are part of the linearization
ivarpar=[ip.periodnullpar,ip.nullparind(:,2)'];
%% extract part of Jacobian
[f,x0]=p_correc_setup(funcs,ptini,ivarpar,method,...
    'remesh_flag',0,'previous',ptini,'output','J',pass_on{:});
ix=dde_ind_from_point(ptini,ivarpar);
[Jext,dum,ieq]=f(x0); %#ok<ASGLU>
vec=@(v)reshape(v,1,[]);
ievar=vec(ix.profile(ip.dim+(1:ip.extdim),:));
ieeq=reshape(ieq.rhsvec,size(ptini.profile,1),[]);
ieeq=vec(ieeq(ip.dim+(1:ip.extdim),:));
J=Jext([ieeq,ieq.extra],[ievar,ix.parameter]);
%% find nullspace
[nullvecs,adj,sv]=dde_svdspaces_lr(J,nulldim,'nullspaces',false); %#ok<ASGLU>
vtemplate=funcs.get_comp(ptini,'vpoint');
vpoint=dde_point_from_x(nullvecs,vtemplate,1:nint);
%% return if only nullspace check is asked for
if strcmp(options.output,'sv')
    varargout={sv,vpoint,ptini};
    return
end
ptini.profile(ievar)=nullvecs(1:length(ievar),1);
ptini.parameter(ip.periodnullpar)=vpoint(1).period(1:length(ip.periodnullpar));
ptini.parameter(ip.nullparind(:,2))=vpoint(1).parameter;
normv2=options.v_scal(ptini,ptini);
ptini.profile(ievar)=ptini.profile(ievar)/sqrt(normv2);
ptini.parameter(ivarpar)=ptini.parameter(ivarpar)/sqrt(normv2);
varargout={ptini,sv,vpoint};
end
%%
function [vtvres,vtvJ]=sys_cond_lincoll_norm(pt,ip,varargin)
%% obtain condition that nullvector/eigenvector has length sqrt(res)
default={'res',1,'period',false};
options=dde_set_options(default,varargin,'pass_on');
vpt=p_axpy(0,pt,[]);
vpt.profile(ip.dim+(1:ip.extdim),:)=pt.profile(ip.dim+(1:ip.extdim),:);
vpt.parameter=pt.parameter;
free_par_ind=[ip.nullparind(:,2)',ip.periodnullpar];
[vtv,vtvJ]=p_dot(vpt,vpt,'free_par_ind',free_par_ind,'period',options.period);
vtvJ=p_axpy(2,vtvJ,[]);
vtvres=vtv-options.res;
end
%%
function pini=dde_collbif_template(ip,pt,omega)
pini=pt;
pini.profile(ip.dim+(1:ip.extdim),:)=0;
pini.parameter([ip.periodnullpar,ip.nullparind(:,2)'])=0;
%% set rotation if type requires it
if isfield(ip,'omega')
    pini.parameter(ip.omega)=omega;
    pini.parameter(ip.period)=pini.period;    
end
end