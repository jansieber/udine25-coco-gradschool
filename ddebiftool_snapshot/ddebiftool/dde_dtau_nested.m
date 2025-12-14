function [dr_out,tauval]=dde_dtau_nested(funcs,x,p,S,dx,dp,dS,tauval)
%% find d(-tau)/d[dx,dp,dT]
%
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