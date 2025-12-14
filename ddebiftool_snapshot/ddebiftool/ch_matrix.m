function y=ch_matrix(funcs,x,p,lam,varargin)
%% Characteristic matrix and its derivatives
% function Delta=ch_matrx(funcs,x,par,lam)
% INPUT:
%   funcs problem functions
%	x steady state solution in R^n (n x 1)
%	p parameter values
%	lam complex number at which characteristic matrix is computed
%   optional named arguments: 'deri', integer (default 0) return derivative
%   of characteristic matrix of order deri wrt lambda
%   'devx', 'devp','order': derivative of ch_matrix wrt x, p deviation
%   (symbols such as {1,'I'} or {2,{2:3}} permitted), of order<=2
%   'dx': if applied to vectors other than identity
%
% function is vectorized such that one can evaluate it for many lambda or
% many x, p, etc.
% OUTPUT: 
%	y characteristic matrix in C^(n x nv)
%
% |x| is assumed to be equilibrium such that delayed values are
% |x(:,ones(1,(ntau+1))| used.
%
%%
default={'dx',{1,'I'},'deri',0,'devx',0,'devp',0,'order',0};
options=dde_set_options(default,varargin,'pass_on');
[lpow_mdx,fac_mdx]=deal(max(1-options.deri,0),double(options.deri<=1));
lpow_dev=options.deri;
args=arg_array_expand([1,2,1],x,p,lam,options.dx,options.devp,0,options.devx,0,0);
args=reshape(args,3,[]);
[x,p,lam,dx]=deal(args{:,1},args{1,2});
vecdim=size(x,2:ndims(x));
[p,lam,dx]=deal(reshape(p,1,size(p,2),[]),...
    reshape(lam,1,[]),reshape(dx,size(dx,1),[]));
[nx,nvec]=size(x);
d=dde_num_delays(funcs)+1; % number of delays incl 0
d_rep=@(arg)repmat(reshape(arg,nx,1,nvec),1,d,1);
x_rep=@(arg)repmat(reshape(arg,1,d,nvec),nx,1,1);
xtau=d_rep(x);
dxtau=d_rep(dx);
taus=dde_taufunvec(funcs,xtau,p,false,true);
vexplt=dxtau.*x_rep(explt_fun(lam,taus,lpow_dev,0)); %nx,d,nvec
ord=options.order;
drhs=@(varargin)funcs.drhs_mf(xtau,p,varargin{:});
lhs_mat=funcs.lhs_matrixfun(nx);
%% char matrix (order 0)
if ord==0
    lam_x=lam(ones(nx,1),:);
    y=fac_mdx*lhs_mat*(dx.*lam_x.^lpow_mdx)-drhs(vexplt,0);
    y=reshape(y,[size(lhs_mat,1),vecdim]);
    return
end
%% derivative of char. matrix wrt x,p (order 1)
[dx2,dp2]=deal(reshape(args{1,3},nx,[]),reshape(args{2,2},1,size(p,2),[]));
dxp2={d_rep(dx2),dp2};
dtaufun=@(order,varargin)funcs.dtau_dir(order,xtau,p,varargin{:});
dtaus=x_rep(dtaufun(1,dxp2{:}));
vdexplt=dxtau.*x_rep(explt_fun(lam,taus,lpow_dev,1));
if ord==1
    y=-drhs(vexplt,0,dxp2{:})-drhs(vdexplt.*dtaus,0);
    y=reshape(y,[size(lhs_mat,1),vecdim]);
    return
end
%% derivative of char. matrix wrt x,p (order 2)
d2taus=x_rep(dtaufun(2,dxp2{:}));
vd2explt=dxtau.*x_rep(explt_fun(lam,taus,lpow_dev,2));
if ord==2
    y=-drhs(vexplt,0,dxp2{:},dxp2{:})-2*drhs(vdexplt.*dtaus,0,dxp2{:})-...
        drhs(vd2explt.*dtaus.^2,0)-drhs(vdexplt.*d2taus,0);
    y=reshape(y,[size(lhs_mat,1),vecdim]);
    return
end
%% higher-order derivative of char. matrix wrt x,p not implemented
error('ch_matrix:order','ch_matrix: derivatives of order %d not implemented',ord);
end
%% y=exp(-z*x)*(-x)^k and its 1st and 2nd derivatives wrt x
% z is 1 x nvec, x is d x nvec, results is d x nvec as z is
% expanded in the first index (z will be lambda, x will be tau)
function y=explt_fun(z,x,k,order)
z=z(ones(size(x,1),1),:);
exz=exp(-z.*x);
if order==0
    y=exz.*((-x).^k);
    return
end
km1=max(k-1,0);
if order==1
    y=exz.*(((-x).^k).*(-z)-k.*((-x).^km1));
    return
end
km2=max(k-2,0);
y=exz.*(((-x).^k).*z.^2+2*z.*k.*((-x).^km1)+k*km1*(-x).^km2);
end