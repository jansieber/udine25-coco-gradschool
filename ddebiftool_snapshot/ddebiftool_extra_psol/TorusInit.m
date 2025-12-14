function varargout=TorusInit(funcs,point,method,varargin)
%% crude initial guess for start of torus bifurcation from Floquet mode
%
%
ip=funcs.ip;
trini=dde_torus_template(funcs,point);
default={'v_scal',@(p,pref)sys_cond_Torus_norm(p,ip.dim,'res',0),'closest',[],...
    'nremove',1,'initmethod','eig','nulldim',2,'align',[],'output','ini'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isempty(options.closest) || strcmp(options.initmethod,'eig')
    [eigval,eigprofile]=mult_crit(funcs.userfuncs,point,method.stability,...
        options.nremove,options.closest);
end
if isempty(options.closest)
    omega=atan2(imag(eigval),real(eigval))/pi;
else
    omega=atan2(imag(options.closest),real(options.closest))/pi;
end
if strcmp(options.output,'sv') ||strcmp(options.initmethod,'svd')|| options.nulldim>2
    out=cell(1,3);
    [out{1:3}]=svd_coll_init(funcs,ip,point,omega,options.nulldim,method.point,...
        'v_scal',options.v_scal,'output',options.output,pass_on{:});
    varargout=out;
    return
end
trini.parameter(ip.omega)=omega;
t=repmat(point.mesh,size(point.profile,1),1);
% convert Floquet multiplier mode to Floquet exponent mode (periodic
% function)
eigprofile=eigprofile.*exp(-log(eigval).*t);
upoint=p_axpy(0,point,[]);
upoint.profile=reshape(real(eigprofile),size(point.profile));
vpoint=upoint;
vpoint.profile=reshape(imag(eigprofile),size(point.profile));
utu=dde_coll_profile_dot(upoint,upoint);
vtv=dde_coll_profile_dot(vpoint,vpoint);
utv=dde_coll_profile_dot(upoint,vpoint);
r=1/sqrt(utu+vtv);
gamma=atan2(2*utv,vtv-utu)/2;
qr=r*(upoint.profile*cos(gamma)-vpoint.profile*sin(gamma));
qi=r*(upoint.profile*sin(gamma)+vpoint.profile*cos(gamma));
trini.profile(ip.null,:)=[qr;qi];
uvpoint=funcs.get_comp(trini,'eigenvector');
sv=NaN(2,1);
varargout={trini,sv,uvpoint};
end
