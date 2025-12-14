function y=sys_rhs_POEV1_var(ip,funcs,order,xx,par,dxx,dpar)
%% use variational problem to compute r.h.s. for fold/BP problem
[nx,nvec]=deal(ip.dim,size(xx,4));
rs=@(v,m)reshape(v,[size(v,1:m),nvec]);
rhs=@(xa,pa)funcs.wrap_rhs(rs(xa(:,:,1,:),2),pa);
rhs_comb=@(varargin)dde_coll_rhs(funcs,varargin{:},'Mxp_incl',false);
beta=rs(par(1,ip.beta,:),1);
T=ones(size(beta));
mxord=2;
x=xx(1:nx,:,1:mxord,:);
v=xx(nx+(1:nx),:,1:mxord,:);
p=par(1,1:ip.nuserpar,:);
q=zeros(size(p));
if isfield(ip,'nullparind')
    q(:,ip.nullparind(:,1),:)=par(:,ip.nullparind(:,2),:);
end
if order==0
    y=[rhs(x,p);...
       rhs_comb(x,p,T, v,q,beta)];
    return
end
dT=zeros(size(beta));
dx=dxx(1:nx,:,1:mxord,:);
dp=dpar(1,1:ip.nuserpar,:);
dq=zeros(size(p));
if isfield(ip,'nullparind')
    dq(:,ip.nullparind(:,1),:)=dpar(:,ip.nullparind(:,2),:);
end
dbeta=rs(dpar(1,ip.beta,:),1);
dv=dxx(nx+(1:nx),:,1:mxord,:);
drhs=@(ord,xa,pa,dxa,dpa)funcs.drhs_dir(ord,rs(xa(:,:,1,:),2),pa,rs(dxa(:,:,1,:),2),dpa);
y1=     drhs(1, x,p,   dx,dp);
y2v=rhs_comb(x,p,T, dv,dq,dbeta);
y2x=rhs_comb(x,p,T,  v, q,beta,  dx,dp,dT);
y=[y1;y2v+y2x];
end
