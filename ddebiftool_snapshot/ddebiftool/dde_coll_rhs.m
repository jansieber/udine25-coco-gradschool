function res=dde_coll_rhs(funcs,varargin)
%% residual for periodic orbit
% arguments after funcs are
%
% xx,p,T,[dx1,dp1,dT1,[dx2,dp2,dT2]],options
%
% xx is nx x d x nderiv x nvec par is 1 x npar x nvec T is 1 x nvec
% followed by deviations of same format dx1, dp1, dT1 (for 1st order
% derivative), dx2 dp2 dT2 (for second order) and options
%
% order 0:
% res=f(x^(k)(t+r)/T^k,p)-M x'(t)/T
%
% r satisfies r=-tau(x(t+r),p)/T, which can be solved recursively in
% ntaucalls (as defined by dde_num_delays) steps, starting from
% r=-tau(x(t+0),p)/T. For constant delays ntaucalls=0. The values of the
% delays are not needed for order==0, since xx contains already the values
% of x at all delayed times.
% 
% order 1: returns directional derivative of res wrt dx1,dp1,dT1, the
% option 'composite' indicates if deviations also apply to delays (value
% true (default)). If the value is false, then delays are treated as
% constant wrt to deviations dx1,dp1,dT1 (resulting in output equal to that
% of funcs.drhs_dir).
%
% order 2: returns mixed-directional derivative of res wrt
% dx1,dp1,dT1,dx2,dp2,dT2. The option 'composite' has now two entries,
% default [true,false], such that the delays are not considered variable
% derivative w.r.t. dx2,dp2,dT2 (this is used for r.h.s functions of
% variational problems for periodic-orbit bifurcation tracking)
%
% output has format nf x nvec
%
%% intermediate variables
% S=1/T,
%
% then the rhs simplifies to res=f(x{k}(t+r)*S^k,p)-M*x'(t)*S
% r  =-S*tau(x{1}(t+r),p)
% 
%% variational equations
%
% 1st order (first index in x,dx is time derivative)
%
% d1S{j}=-1/T^2 dTj, d1(S^k){j}=-k/T^(k+1) dTj (j=1,2)
% d2S=2/T^3 dT1 dT2
%
% dr{j,1}  =-d1S{j}*tau-S*dtau(x(t+r),p)[dxj(t+r),       dpj]
% dr{j,i+1}= dr{j,1}-   S*dtau(x(t+r),p)[x'(t+r)*dr{j,i}, 0 ]
% for i=2..ntaucalls (nesting number of delays)
%
% 0 = M x'(t)dTj/T^2 - M dxj'(t)/T + df(y,p)[dy,dp]
% where y=x^(k)(t+r)*S^k
% dy=dxj^(k)(t+r)*S^k + x^(k+1)(t+r) dr{j,end}*S^k +x^(k)(t+r)*d(S^k){j} 
% k=0..mxd, r has d components (first is 0), j=1,2. t in (0,1)
%
%% permit arguments as cell or as list
default={'Mxp_incl',true,'tau',NaN,'composite',[true,false]};
[order,x,p,T,dx1,dp1,dT1,dx2,dp2,dT2,opts]=select_args(default,varargin{:});
[nx,d,nderivp1]=size(x,1:3);
[x,p,T,vecdim,nvec]=arg_flatten([3,2,1],x,p,T);
M=funcs.lhs_matrixfun(nx);
rhsfun=@(varargin)funcs.drhs_dir('incl_deriv',varargin{:});
ntaum1=dde_num_delays(funcs);
assert(d==ntaum1+1); % check if delays are ok
rep=@(v,shape,rep)repmat(reshape(v,shape{:}),rep);
srepxd=@(v)repmat(reshape(v,1,1,nderivp1,nvec),[nx,d,1,1]);
[S,S_k] = S_T(T,nderivp1,nvec);
S_k_xm=srepxd(S_k); % S^k (k=0..mxd+1)
%%
D= @(x,deg)reshape(x(:,1,1+deg,:),nx,nvec);
Ddm=@(x,deg)reshape(x(:,:,1+deg,:),nx,d,length(deg),nvec);
%% residual
if order==0
    lhs=0;
    if opts.Mxp_incl
        lhs=M*(D(x,1).*repmat(S,nx,1));
    end
    res=-lhs+rhsfun(0,x(:,:,1:nderivp1,:).*S_k_xm,p);
    res=reshape(res,[size(res,1),vecdim]);
    return
end
%% variational equations
is_composite=opts.composite;
%% expand/reshape arrays
[dx1,dp1,dT1]=arg_flatten([3,2,1],dx1,dp1,dT1);
[d1S1,d1S1_k]=dS1_T(T,dT1,nderivp1,nvec);
d1S1_k_xm=srepxd(d1S1_k);
dp0=zeros(size(dp1));
rs=@(v,n)reshape(v,[size(v,1:n),nvec]);
%% determine d[x^(k)(t+r)*S^k]
Ddp=@(ua,i)Ddm(ua,i+(0:nderivp1-2));
Edp=@(ua)Ddp(ua,0);
dx_r{1}=Edp(x.*d1S1_k_xm+dx1.*S_k_xm);
%% determine dr recursively
if is_composite(1)
    [dr{1},tauval]=dde_dtau_nested(funcs,x,p,S,dx1,dp1,d1S1,opts.tau);
    d1rk1x=rep(dr{1},{1,d,1,nvec},[nx,1,nderivp1-1,1]);
    dx_r{1}=dx_r{1}+Ddp(x,1).*Edp(d1rk1x).*Edp(S_k_xm);
end
%% combine rhs and lhs
if order==1
    lhs=M*(D(x,1).*rs(d1S1_k_xm(:,1,2,:),1));
    if opts.Mxp_incl % dx1' is only needed here
        lhs=lhs+M*rs(dx1(:,1,2,:).*S_k_xm(:,1,2,:),1);
    end
    df=rhsfun(1,Edp(x.*S_k_xm),p,dx_r{1},dp1);
    res=df-lhs;
    res=reshape(res,[size(res,1),vecdim]);
    return
end
%% 2nd derivative of S^k
[d2S,d2S_k]  =dS2_T(T,dT1,dT2,nderivp1,nvec);
%% 1st derivative of S^k wrt 2nd deviation
[d1S2,d1S2_k]=dS1_T(T,    dT2,nderivp1,nvec);
d1S2_k_xm=srepxd(d1S2_k); % d[S^k]/dT [dT2]
d2S_k_xm=srepxd(d2S_k); % d^2[S^k][dT1][dT2]
dx_r{2}=Edp(x.*d1S2_k_xm+dx2.*S_k_xm);
%% find first derivative of delays wrt to 2nd deviation
if is_composite(2)
    [dr{2},tauval]=dde_dtau_nested(funcs,x,p,S,dx2,dp2,d1S2,tauval);
    d1rk2x=rep(dr{2},{1,d,1,nvec},[nx,1,nderivp1-1,1]);
    dx_r{2}=dx_r{2}+Ddp(x,1).*Edp(d1rk2x).*Edp(S_k_xm);
end
%% determine d2r recursively
if any(is_composite)
    d2r=dr2_by_dxpT(funcs,is_composite,...
        x,p,S,dx1,dp1,d1S1,dx2,dp2,d1S2,d2S,tauval,dr);
    d2rkx=rep(d2r,{1,d,1,nvec},[nx,1,nderivp1-1,1]);
end
%% determine d^2[x^(k)(t+r)*/T^k]
mxd=1+double(any(is_composite))+double(all(is_composite));
Ddp=@(ua,i)Ddm(ua,i+(0:nderivp1-mxd));
Edp=@(ua)Ddp(ua,0);
d2x_r=...
      Edp( x   ).*Edp(d2S_k_xm )+...
      Edp(dx1  ).*Edp(d1S2_k_xm)+...
      Edp(dx2  ).*Edp(d1S1_k_xm);
if is_composite(1)
    d2x_r=d2x_r+...
        Ddp(x,1).*Edp(d1S2_k_xm).*Edp(d1rk1x)+...
        Ddp(dx2,1).*Edp( S_k_xm).*Edp(d1rk1x);
end
if is_composite(2)
    d2x_r=d2x_r+...
        Ddp(x,1).*Edp(d1S1_k_xm).*Edp(d1rk2x)+...
        Ddp(dx1,1).*Edp( S_k_xm).*Edp(d1rk2x);
end
if all(is_composite)
    d2x_r=d2x_r+Ddp(x,2).*Edp(S_k_xm).*Edp(d1rk1x.*d1rk2x);
end
if any(is_composite)
    d2x_r=d2x_r+Ddp(x,1).*Edp(S_k_xm).*Edp(d2rkx);
end
%% combine rhs and lhs
rhsmf=@(varargin)funcs.drhs_mf('incl_deriv',Edp(x.*S_k_xm),p,varargin{:});
if order==2
    df2 =rhsmf(dx_r{1},dp1, dx_r{2},dp2);
    df11=rhsfun(1,Edp(x.*S_k_xm),p,d2x_r,dp0);
    lhs=dx1(:,1,2,:).*d1S2_k_xm(:,1,2,:)+...
        dx2(:,1,2,:).*d1S1_k_xm(:,1,2,:)+...
         x( :,1,2,:).*d2S_k_xm( :,1,2,:);
    res=-M*rs(lhs,1)+df2+df11;
    res=reshape(res,[size(res,1),vecdim]);
    return
end
end
%%
function varargout=select_args(default,varargin)
if iscell(varargin{1})
    args=varargin{1};
    optscell=varargin(2:end);
else
    args=varargin;
    optstart=find(cellfun(@ischar,args),1,'first');
    if isempty(optstart)
        optscell={};
    else
        optscell=varargin(optstart:end);
        args=args(1:optstart-1);
    end
end
args=args(:).';
order=length(args)/3-1;
args=[args,cell(1,9-length(args))];
opts=struct(default{:});
for i=1:2:length(optscell)
    if isfield(opts,optscell{i})
        opts.(optscell{i})=optscell{i+1};
    end
end
varargout=[{order},args,{opts}];
end
%% S=1/T, and powers
function [S,S_k]=S_T(T,nderivp1,nvec)
S=1./T;
kpow=0:nderivp1-1;
kpowrep=repmat(kpow,[1,1,nvec]);
S_k=repmat(reshape(S,1,1,nvec),[1,nderivp1,1]).^kpowrep;
end
%% 1st-order derivative of S^k wrt T, deviation dT
function [dS,dS_k]=dS1_T(T,dT,nderivp1,nvec)
kpow=repmat(0:nderivp1-1,[1,1,nvec]);
dS=(-1./T.^2).*dT;
rep=@(v)repmat(reshape(v,1,1,nvec),[1,nderivp1,1]);
dS_k=-kpow./rep(T).^(kpow+1).*rep(dT);
end
%% 2nd-order derivative of S^k wrt T, deviations dT1,dT2
function [d2S, d2S_k] = dS2_T(T,dT1,dT2,nderivp1,nvec)
d2S=(2./T.^3).*dT1.*dT2;
kpow=repmat(0:nderivp1-1,[1,1,nvec]);
srepp=@(v)repmat(reshape(v,1,1,nvec),[1,nderivp1,1]);
d2S_k=kpow.*(kpow+1)./srepp(T).^(kpow+2).*srepp(dT1.*dT2);
end
%% find dr (1st order)
% x is nx x d x nderivp1 x nvec, p is 1 x npar x nvec, S is 1 x nvec
% deviations have same format
%
% argument S=1/T, dS=[dS/dT]dT, input tauval may contain precomuted delays
% tau (>0), if not they are returned.
%
%% 
% outputs are derivatives of -tau/T
%
% r=-S*taufun(x(t+r),p), implicitly differentiated wrt. x,p,S and
% recursively solved with loop of length ntaucalls
function [dr_out,tauval]=dde_dtau_nested(funcs,x,p,S,dx,dp,dS,tauval)
[nx,d,nderivp1,nvec]=size(x); %#ok<ASGLU>
Ed=@(x)reshape(x(:,1:d,1,:),nx,d,nvec);
Dd=@(x,deg)reshape(x(:,:,1+deg,:),nx,d,nvec);
srepd=@(sa)repmat(sa,d,1);
drepx=@(dra)repmat(reshape(dra,1,d,nvec),[nx,1,1]);
ntaucalls=dde_num_delays(funcs,'ntaucalls');
taufun=@(order,varargin)funcs.dtau_dir(order,Ed(x),p,varargin{:});
if isnan(tauval)
    tauval=taufun(0);
end
dr=NaN(d,nvec,ntaucalls);
[S,dS]=deal(srepd(S),srepd(dS));
dr(:,:,1)=-dS.*tauval-S.*taufun(1,Ed(dx),dp);
dp0=zeros(size(dp));
for k=2:ntaucalls
    dr(:,:,k)=dr(:,:,1)-S.*taufun(1,Dd(x,1).*drepx(dr(:,:,k-1)),dp0);
end
dr_out=dr(:,:,end);
end
%%
function d2r_out=dr2_by_dxpT(funcs,is_composite,x,p,S,dx1,dp1,d1S1,dx2,dp2,d1S2,d2S,...
    tauval,dr)
[nx,d,nderivp1,nvec]=size(x); %#ok<ASGLU>
npar=size(dp1,2);
Ed=@(x)reshape(x(:,1:d,1,:),nx,d,nvec);
Dd=@(x,deg)reshape(x(:,:,1+deg,:),nx,d,nvec);
rep=@(v,shape,rep)repmat(reshape(v,shape{:}),rep);
drepx=@(dra)rep(dra,{1,d,nvec},[nx,1,1]);
ntaucalls=dde_num_delays(funcs,'ntaucalls');
taufun=@(order,varargin)funcs.dtau_dir(order,Ed(x),p,varargin{:});
taumf=@(varargin)funcs.dtau_mf(Ed(x),p,varargin{:});
if isnan(tauval)
    tauval=taufun(0);
end
srepd=@(s)repmat(s,d,1);
sreppar=@(s)rep(s,{1,1,nvec},[1,npar,1]);
srepxd0=@(s)rep(s,{1,1,nvec},[nx,d,1]);
[Sxd,d1S1xd,d1S2xd]=deal(srepxd0(S),srepxd0(d1S1),srepxd0(d1S2));
dp0=zeros(size(dp1));
tau1_xarg=d1S1xd.*Ed(dx2)+d1S2xd.*Ed(dx1);
tau2_xargs={Ed(dx1),Ed(dx2)};
tau1_parg=sreppar(d1S1).*dp2+sreppar(d1S2).*dp1;
if is_composite(1)
    tau1_xarg=tau1_xarg+(d1S2xd.*Dd(x,1)+Sxd.*Dd(dx2,1)).*drepx(dr{1});
    tau2_xargs{1}=tau2_xargs{1}+Dd(x,1).*drepx(dr{1});
end
if is_composite(2)
    tau1_xarg=tau1_xarg+(d1S1xd.*Dd(x,1)+Sxd.*Dd(dx1,1)).*drepx(dr{2});
    tau2_xargs{2}=tau2_xargs{2}+Dd(x,1).*drepx(dr{2});
end
if all(is_composite)
    tau1_xarg=tau1_xarg+Sxd.*Dd(x,2).*drepx(dr{1}).*drepx(dr{2});
end
d2r=NaN(d,nvec,ntaucalls);
d2r(:,:,1)=-srepd(d2S).*tauval-taufun(1,tau1_xarg,tau1_parg)-...
    srepd(S).*taumf(tau2_xargs{1},dp1,tau2_xargs{2},dp2);
for k=2:ntaucalls
    d2r(:,:,k)=d2r(:,:,1)-srepd(S).*taufun(1,Dd(x,1).*drepx(d2r(:,:,k-1)),dp0);
end
d2r_out=d2r(:,:,end);
end
