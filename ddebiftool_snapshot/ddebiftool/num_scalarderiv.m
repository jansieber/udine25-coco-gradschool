%% Single-variable derivative using numerical finite difference
%
% $$ y=\frac{d^m}{dx^m}f(x)\vert_{x=0} $$
%
% arguments: |rhs=f|, |m=dord|. Additonal
% arguments for numerical differentiation:
% |fac| (1.99) ratio between distances of interpolation points from 0:
% interpolation points are:
% -fac^(n-1)*h,...-fac*h,-h,0,h,h*fac,...,h*fac^(n-1),
% |h| (default eps^(1/(m+1))) scaling in spacing and output matrix
%
%% Input
%
% * |rhs_inp| r.h.s, y=f(x), which gets differentiated
% * |dord| order of requested derivative
%
%% Optional inputs
%
% * |h| (default |eps^(1/(dord+1))|) initial scaling in spacing and output matrix
% * |fac| (1.99) ratio between successive deviations
% * |isvectorized| can rhs be called with several x?
% * |maxit| (20) maximal number of adjustments of h by fac
% * |output| (1): switch to 2 to swap order of first two outputs  
%
%% Output
%
% * |y| derivative of order |dord|
% * |y2| derivative of order |dord|, lower order estimate
% * |h| spacing/scaling used in final estimate
%%
function [y,y2,h]=num_scalarderiv(rhs_inp,dord,varargin)
% process optional arguments
default={'h',eps^(1/(dord*2)),'fac',1.99,'isvectorized',true,'maxit',10,'output',1};
options=loc_set_opts(default,varargin);
%% approximation order is 2,5,6,9,10,... for derivatives of order 1,2,3,4,5,..
fac=options.fac;
deriv=diff_bary_wt(dord,fac,options.h,dord);
deriv.x=deriv.x*options.h;
% wrapper permitting to assume vectorization
rhs=@(x)rhs_vec(rhs_inp,x,options.isvectorized); 
% apply linear combination to inerpolants to obtain approx derivative
D=@(ind,deriv,h)(deriv.f0*deriv.D0+deriv.fvals(:,2*ind+(1:length(deriv.D)))*deriv.D)/h^dord;
%% Compute initially 3 approximations for derivative for h/fac, h h*fac
% We first check which two approximations are closer to each other.
% Depending on that we decrease (or increase) h interatively until the
% closeness relation changes.
deriv.f0=rhs(0);
szfv=size(deriv.f0);
deriv.x=[...
    deriv.x(:,1)/fac, deriv.x, deriv.x(:,end)*fac];
deriv.fvals=rhs(deriv.x);
h=options.h;
df=[D(0,deriv,h/fac),D(1,deriv,h),D(2,deriv,h*fac)];
err=max(abs(diff(df,[],2)),[],1);
err_dir=err(1)<err(2);
for it=1:options.maxit
    if err_dir ~=(err(1)<err(2))
        break
    end
    if err(1)<err(2)
        deriv.x=[deriv.x(:,1)/options.fac,deriv.x(:,1:end-1)];
        deriv.fvals=[rhs(deriv.x(:,1)),deriv.fvals(:,1:end-2)];
        h=h/fac;
    elseif err(1)>=err(2)
        deriv.x=[deriv.x(:,2:end),deriv.x(:,end)*options.fac];
        deriv.fvals=[deriv.fvals(:,3:end),rhs(deriv.x(:,end))];
        h=h*fac;
    end
    df=[D(0,deriv,h/fac),D(1,deriv,h),D(2,deriv,h*fac)];
    err=max(abs(diff(df,[],2)),[],1);
end
%% Final output
% apply Richardson extrapolation to the approximations h,h*fac
% "lower order" estimate is either the approximation for h/fac or h*fac,
% depending on which one is further away from the final output.
rfac=1/fac^deriv.approx_order;
extrapmat=[[1;-rfac;0], [0;1;-rfac]]/(1-rfac);
df_extrap=df*extrapmat;
y=reshape(df_extrap(:,2),szfv);
err1=norm(df_extrap(:,2)-df(:,1),inf);
err3=norm(df_extrap(:,2)-df(:,3),inf);
if err1>err3
    y2=reshape(df(:,1),szfv);
else
    y2=reshape(df(:,3),szfv);
end
if options.output==2
    [y,y2]=deal(y2,y);
end
end
%% Matrix corresponding to kth derivative of interpolant around 0
%%
function ret=diff_bary_wt(order,fac,h,npoints)
%
%% Inputs:
% * order: integer, order k of derivative
% * fac: spacing of interpolation points is 
%       -fac^(n-1)*h,...-fac*h,-h,0,h,h*fac,...,h*fac^(n-1)
% * h: scaling in spacing and output matrix
% * npoints: number of interpolation points to each side
%
%% Output structure ret. Fields:
% * D (1 x n) coefficients for non-zero interpolcation points
% * x (2 x n) interpolcation points to each side: x(:,k+1)=h*fac^k*[-1;1]
% * D0: coeficient at 0
% * deriv: order of derivative
% * approx_order: interpolcation order of formula (2*order for odd,
% 2*order+1 for even order derivatives)
%
% $Id: nmfm_Dbary.m 309 2018-10-28 19:02:42Z jansieber $
%
%% interpolcation points
facpow=fac.^(0:npoints-1);
facpowsgn=[-facpow;facpow];
x=[0,facpowsgn(:)']*h;
nx=length(x);
%% create barycentric weights wi
[xi1,xi2]=ndgrid(x,x);
dxi=xi1-xi2+eye(nx);
w=1./prod(dxi,2)';
on=ones(1,nx);
wrep=w(on,:);
xreph=x(on,:);
%% Barycentric formula for 1st derivative at interpolation points
% gives matrix D (see paper by Berrut, Trefethen ,SIAM Review 46(3))
denom=(xreph-xreph');
denom(1:nx+1:end)=Inf;
D=-wrep./wrep'./denom;
Drsum=sum(D,2);
D(1:nx+1:end)=-Drsum;
%% D^order is (higher) order derivative interpolation
D=D^order;
D=D(1,:)';
ret=struct('D',D(2:end),'x',reshape(x(2:end),2,[]),'D0',D(1),...
    'deriv',order,'approx_order',order*2+1);
if mod(order,2)
    ret.D0=0;
    ret.approx_order=ret.approx_order-1;
end
end
%% wrapper for vectorized call to rhs
function f=rhs_vec(rhs,x,isvec)
if isvec
    f=rhs(x(:)');
else
    f=cell2mat(arrayfun(@(i){rhs(x(i))},1:numel(x)));
end
f=reshape(f,[],numel(x));
end
%%
function opts=loc_set_opts(default,args)
opts=struct(default{:});
args=reshape(args,2,[]);
for i=1:size(args,2)
    if isfield(opts,args{1,i})
        opts.(args{1,i})=args{2,i};
    end
end
end
