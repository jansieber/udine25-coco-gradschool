function res=dde_psol_sysvar1(xx,p,dxx,dp,dTs,funcs,tau,ncalls)
%% variational problem for DDE, scaled to psol interval [0,1]
[nx,ntau,nderivs,nvec]=size(xx); %#ok<ASGLU>
%% name x, x', dx
rs=@(x)reshape(x,nx,ntau,[]);
[E0x,E1x,E0dx]=deal(rs(xx(:,:,1,:)),rs(xx(:,:,2,:)),rs(dxx(:,:,1,:)));
%% find delays 
if nargin<8
    [tau,ncalls]=dde_taufunvec(funcs,E0x,p);
end
%% repeat parameters and deviations if necessary
dtau=@(dxa,dpa)funcs.dtau_dir(1,E0x,p,E0dx,dpa); % first order, include zero delay
rrep=@(v,shape,rep)repmat(reshape(v,shape{:}),rep);
M=funcs.lhs_matrixfun(nx);
diag0m=@(xv,tv)xv.*rrep(tv,{1,ntau,nvec},[nx,1,1]);
dtaudp=dtau(zeros(size(dxx(:,1:ntau,:))),dp);
dTsmat=rrep(dTs,{1,nvec},[ntau,1]);
% dy=dx+x'*tau*dT-x'*(dtau/dp)dp for constant delays
dy0ini=E0dx+diag0m(E1x,tau.*dTsmat)-diag0m(E1x,dtaudp); 
%% for state-dependent delays
dy=dy0ini;
for i=1:ncalls
    dy0ini=-diag0m(E1x,dtau(dy0ini,zeros(size(dp))));
    dy=dy+dy0ini;
end
%% linearization of r.h.s.
xp=reshape(E1x(:,1,:),nx,nvec);
dTsvec=rrep(dTs,{1,nvec},[nx,1]);
res=M*(xp.*dTsvec)+funcs.drhs_dir(1,xx(:,1:ntau,:),p,dy,dp);
end
%% add higher-order derivative args (only non-empty if nderivs>1)
% rrrep=@(v,shp1,rep,shp2)reshape(repmat(reshape(v,shp1{:}),rep),shp2{:});
% mxf=nderivs-1;
% rgk=ntau+1:ntau_orig;
% [Ek1x,Ek0dx]=deal(xx(:,ntau+rgk,:),dxx(:,rgk,:));
% diagkm=@(xv,tv)xv.*rrep(tv,{1,ntau,nvec},[nx,mxf,1]);
% deriv_vec=rrrep((1:mxf),{1,1,mxf,1},[nx,ntau,1,nvec],{nx,ntau*mxf,nvec});
% dyk=Ek0dx+diagkm(Ek1x,tau.*dTsmat)-diagkm(Ek1x.*deriv_vec,dTsmat)-diagkm(Ek1x,dtau(dy0,dp));
%
% dy=reshape(cat(2,dy0,dyk),nx,ntau_orig,nvec);
