function y=rot_rhs(A,expA,fcn,lhs_matrix,ord,xx,p,dx,dp)
%% right-hand side in rotating coordinates
%
% $Id: rot_rhs.m 369 2019-08-27 00:07:02Z jansieber $
%
omega=reshape(p(1,end,:),1,[]);
userpar=p(1,1:end-1,:);
[nx,d,nvec]=size(xx);
rs=@(x)reshape(x(:,1,:),nx,nvec);
tau=fcn.dtau_dir(0,xx,userpar); % delays including 0
expAmat=Momx(expA,omega,'matrix',-tau);
xxrot=reshape(expAmat*xx(:),nx,d,nvec);
if ord==0
    y0=fcn.wrap_rhs(xxrot,userpar);
    y1=Momx(lhs_matrix*A,omega,rs(xx));
    y=y0-y1;
    return
end
%% 1st derivative
[dom,dpar]=deal(reshape(dp(1,end,:),1,[]),dp(1,1:end-1,:));
dtau=fcn.dtau_dir(1,xx,userpar,dx,dpar);
dxrot=reshape(expAmat*dx(:),nx,d,nvec);
domtau_mat=Momx(A,dom,'matrix',-tau)+Momx(A,omega,'matrix',-dtau);
dexpAxrot=dxrot+reshape(domtau_mat*xxrot(:),nx,d,nvec);
if ord==1
    y0=fcn.drhs_dir(1,xxrot,userpar,dexpAxrot,dpar);
    y1=Momx(lhs_matrix*A,dom,rs(xx))+Momx(lhs_matrix*A,omega,rs(dx));
    y=y0-y1;
    return
end
end
%% multiply M*omega*tau or expM(omega*tau) with x
% if ix is not numeric, return matrix is sparse_blkdiag form
function y=Momx(M,om,x,tau)
if nargin>3 && ~isempty(tau)
    omtau=om(ones(size(tau,1),1),:).*tau;
else
    omtau=om;
end
if ~isnumeric(M)
    Msp=sparse_blkdiag(M(reshape(omtau,1,[])));
else
    omtau=repmat(reshape(omtau,1,1,[]),[size(M),1]);
    Msp=sparse_blkdiag(repmat(M,[1,1,size(omtau,3)]).*omtau);
end
if isnumeric(x)
    y=reshape(Msp*x(:),size(M,1),[]);
else
    y=Msp;
end
end
