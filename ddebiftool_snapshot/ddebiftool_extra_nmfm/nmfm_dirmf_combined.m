function [y,cfpol]=nmfm_dirmf_combined(funcs,x0,p0,order,taus,xdevs,pdevs,varargin)
%% compute derivative of order for composite functional in single direction
% at equilibrium
%
% [d^k/dh^k] F(x+hv,p+hq)=[d^k/dh^k] f(y(h),p+hq)
%
% assume y=sum yj h^j/j!, j=0..k, determine y first:
%
% yj=j [d^(j-1)/dh^(j-1)] v(-tau(y(h),p+hq))
%
% y0=x(-tau)=pt.x (equilibrium)
% y1=v(-tau)
%
% $Id: nmfm_dirmf_combined.m 309 2018-10-28 19:02:42Z jansieber $
%%
default={'print',0,'cfpol',[]};
options=dde_set_options(default,varargin,'pass_on');
% coeffs for chain rule for up to maxorder
if isempty(options.cfpol)
    options.cf=dde_chainrule_combinatorics(order);
    for k=length(options.cf):-1:1
        cfpol{k}=dde_pol_from_chain(options.cf(k));
    end
else
    cfpol=options.cfpol;
end
ndim=size(x0,1);
nvec=size(x0,3);
npvec=size(p0,3);
assert(size(p0,1)==1 && (npvec==1 || npvec==nvec));
%% compute delays at equilibrium
ntaux=length(taus);
np=size(p0,2);
% determine deviations v of solutions and their derivatives in time, and
% deviation q of parameters. v(:,j,k,l) is (l-1)st deriv of v_j at tau_k
% (starting to count tau at tau_1=0).
nz=ndim*ntaux;
z=zeros(nz,nvec,order+1);
z(1:nz,:,1)=reshape(x0(:,:,:,1),nz,nvec);
z(nz+(1:np),:,1)=reshape(p0,np,nvec);
z(1:nz,:,2)=reshape(xdevs(:,:,:,1),nz,nvec);
z(nz+(1:np),:,2)=reshape(pdevs,np,nvec);
taus=repmat(taus,1,nvec);
dtauvals=cat(3,reshape(taus,ntaux,nvec,1),zeros(ntaux,nvec,order-1));
f_dtau=@(order,y0,dy)dtaufunc(order,y0,dy,funcs,ndim,ntaux,nvec,npvec);
for i=2:order
    %% compute d^i/dh^i [tau(y(h),p+hq)]
    dtauvals(:,:,i)=nmfm_chainrule(f_dtau,z(:,:,1:i),i-1,'cf',cfpol{i-1});
    %% compute d^(i+1)/dh^(i+1) y(h)
    dvdti=@(order,tau0,dtau)dvfunc(order,dtau,xdevs);
    z(1:nz,:,i+1)=nmfm_chainrule(dvdti,dtauvals(:,:,1:i),i-1,'cf',cfpol{i-1});
    z(1:nz,:,i+1)=i*z(1:nz,:,i+1);
end
f_drhs=@(order,y0,dy)drhsfunc(order,y0,dy,funcs.drhs_dir,ndim,ntaux,nvec,npvec);
y=nmfm_chainrule(f_drhs,z,order,'cf',cfpol{order});
end
%% derivative(s) of v(-tau)
function zval=dvfunc(order,dtau,v)
ndim=size(v,1);
sgn=mod(order-1,2)*2-1;
dtau=reshape(dtau,[1,size(dtau)]);
zval=sgn*v(:,:,:,order+1).*(dtau(ones(1,ndim),:,:)).^order;
zval=reshape(zval,[],size(v,3));
end
%% derivative(s) of sys_tau
function dtauval=dtaufunc(order,y0,dy,funcs,ndim,ntaux,nvec,nvecp)
x0=reshape(y0(1:ndim*ntaux,:),ndim,ntaux,nvec);
p0=reshape(y0(ndim*ntaux+1:end,1:nvecp),1,[],nvecp);
dx =reshape(dy(1:ndim*ntaux,:),size(x0));
dp =reshape(dy(ndim*ntaux+1:end,1:nvecp),size(p0));
dtauval=funcs.dtau_dir(order,x0,p0,dx,dp);
[ntaum1,max_rhs_xderiv]=dde_num_delays(funcs); %#ok<ASGLU>
dtauval=repmat(dtauval,max_rhs_xderiv+1,1);
end
%% derivative(s) of sys_rhs
function drhsval=drhsfunc(order,y0,dy,rhsdirderi,ndim,ntaux,nvec,nvecp)
x0=reshape(y0(1:ndim*ntaux,:),ndim,ntaux,nvec);
p0=reshape(y0(ndim*ntaux+1:end,:),1,[],nvecp);
v =reshape(dy(1:ndim*ntaux,:),size(x0));
q =reshape(dy(ndim*ntaux+1:end,:),size(p0));
drhsval=rhsdirderi(order,x0,p0,v,q);
end
