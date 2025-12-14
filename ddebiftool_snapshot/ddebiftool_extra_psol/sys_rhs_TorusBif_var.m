function y=sys_rhs_TorusBif_var(ip,funcs,order,xx,par,dxx,dpar)
%% rhs for torus bifurcation of periodic orbits in (SD-)DDEs
[nx,ntau,nvec]=deal(ip.dim,ip.orig_ntau+1,size(xx,4));
rs=@(v,m)reshape(v,[size(v,1:m),nvec]);
rep=@(v,shape,rep)repmat(reshape(v,shape{:}),rep);
drepx1=@(v)rep(v,{1,ntau,1,nvec},[nx,1,1,1]);
%% determine rotation
[omega,Period]=deal(par(1,ip.omega,:),par(1,ip.period,:));
rho=pi*rs(omega./Period,1);
%% extract user states and parameters
x=xx(ip.xrg,:,1:2,:);
x0=rs(x(:,:,1,:),2);
p=par(1,1:ip.nuserpar,:);
%% variational problem does not include variations wrt period or parameters
[dp0,dT0,T1]=deal(zeros(size(p)),zeros(1,nvec),ones(1,nvec));
%% delay values enter deviations explicitly, so need to be recomputed
tau=dde_taufunvec(funcs,x0,p);
%% rotation of variation
e_irhotau=drepx1(exp(-1i*rho(ones(ntau,1),:).*tau)); % ntau x nvec
z=xx(ip.re,:,1,:)+1i*xx(ip.im,:,1,:);
%% z(t-tau) exp(-i rho tau)
z_ertau=z.*e_irhotau; 
M=funcs.lhs_matrixfun(nx);
irz=1i*rho(ones(nx,1),:).*rs(z(:,1,1,:),1);
rhs_nested=@(varargin)dde_coll_rhs(funcs,varargin{1:end-1},...
    'Mxp_incl',false,'tau',varargin{end});
rhs_nested_c=@(varargin)complex_eval(rhs_nested,4,x,p,T1,varargin{:},tau);
if order==0
    df=rhs_nested_c(z_ertau,dp0,dT0)-M*irz;
    y=[funcs.wrap_rhs(x0,p);...
        real(df);imag(df)];
    return
end
%% deviations 
[dp2,dom,dPeriod]=deal(dpar(1,1:ip.nuserpar,:),dpar(1,ip.omega,:),dpar(1,ip.period,:));
drho=pi*rs(dom./Period,1)-pi*rs(omega./Period.^2.*dPeriod,1);
dx=dxx(ip.xrg,:,1:2,:);
dz0=dxx(ip.re,:,1,:)+1i*dxx(ip.im,:,1,:);
dx20=rs(dx(:,:,1,:),2);
%% 
irdz=1i*rho( ones(nx,1),:).*rs(dz0(:,1,1,:),1);
idrz=1i*drho(ones(nx,1),:).*rs(  z(:,1,1,:),1);
%% find derivatives of delays wrt dx2, dp2
rrepxd=@(r)rep(r,{1,1,1,nvec},[nx,ntau,1,1]);
dtaux=drepx1(funcs.dtau_dir(1,x0,p,dx20,dp2));
[rhox,drhox]=deal(rrepxd(rho),rrepxd(drho));
df1xarg=(dz0-1i*z.*drepx1(tau).*drhox-1i*rhox.*z.*dtaux).*e_irhotau;
z_ertau(:,:,2,:)=0; % time derivative of z_ertau should not be needed!
y1=funcs.drhs_dir(1,x0,p,dx20,dp2);
y2d1=rhs_nested_c(df1xarg,dp0,dT0)-M*irdz-M*idrz;
y2d2=rhs_nested_c(z_ertau,dp0,dT0,  dx,dp2,dT0);
if order==1
    y=[y1;real(y2d1+y2d2);imag(y2d1+y2d2)];
end
end
%% evaluate fcn on complex deviations
function res=complex_eval(fcn,complexargs,varargin)
args=varargin;
iscomplex=false(1,length(args));
iscomplex(complexargs)=true;
for i=1:length(args)
    if iscomplex(i)
        args{i}=cat(ndims(args{i}),real(args{i}),imag(args{i}));
    else
        args{i}=cat(ndims(args{i}),args{i},args{i});
    end
end
res2=fcn(args{:});
nvec=size(res2,2)/2;
res=res2(:,1:nvec)+1i*res2(:,nvec+(1:nvec));
end

