function y=int_delay_chain_element(gfun,msh,order,varargin)
%% approximate algebraic equation for integral distributed delay
%
% y(t-r,1)=x0(t)+int_0^r g(s,x(t-s),p) dr
% y(t-r,2)=g(r,x(t-r),p)  for r=0..tau
%
% discretized
%
% y(:,i,:,1)=intwJ*[x0(:,t);tau*g(t(2:end)*tau,xx(:,t-t(2:end)*tau),p)]
%
% where intwJ and t(i) are predefined integration weights and nodes on the
% interval [0,1], if the interval is finite, or on [0,infty] if the
% interval is infinite (Laguerre or generalised Laguerre approximation)
%
% derivatives up to order 2 provided
%%
assert(all(order<=2));
nbaseargs=4;
[xx,p,tau,x0]=deal(varargin{1:nbaseargs});
[nw,nx,nvec]=deal(length(msh.t),size(xx,1),size(xx,3));
ng=size(msh.wJ,2)/nw;
tau=reshape(tau,1,nvec);
bdrep=tau(ones(ng*nw,1),:);
s= reshape(reshape(msh.t,[],1)*tau,1,1,nw*nvec);
xs=reshape(xx,nx,1,nw*nvec);
ps=reshape(repmat(p,1,nw,1),1,size(p,2),nw*nvec);
g=@(ord,varargin)gfun(ord,s,xs,ps,varargin{:});
gval=reshape(g(0),ng*nw,nvec);
[y,done]=deal({},false(size(order)));
[y,done]=combine_g(x0,bdrep.*gval,gval,0,order,y,done,msh);
if all(done)
    return
end
[dxx,dp,dtau,dx0]=deal(varargin{nbaseargs+(1:nbaseargs)});
dtau=reshape(dtau,1,nvec);
dbdrep=dtau(ones(ng*nw,1),:);
ds= reshape(msh.t(:)*dtau,1,1,nw*nvec);
dxs=reshape(dxx,nx,1,nw*nvec)  ;
dps=reshape(repmat(dp,1,nw,1),1,size(dp,2),nw*nvec);
d1gval1=reshape(g(1,ds,dxs,dps),ng*nw,nvec);
[y,done]=combine_g(dx0,dbdrep.*gval+bdrep.*d1gval1,d1gval1,1,order,y,done,msh);
if all(done)
    return
end
dx2s=reshape(varargin{2*nbaseargs+1},nx,1,nw*nvec);
[ds0,dps0,dx00]=deal(zeros(size(ds)),zeros(size(dps)),zeros(size(dx0)));
d1gval2=reshape(g(1,ds0,dx2s,dps0),ng*nw,nvec);
d2gval= reshape(g(2,ds, dxs, dps ),ng*nw,nvec);
dg2=d1gval2+d2gval;
dg2total=2*dbdrep.*d1gval1+bdrep.*dg2;
[y,done]=combine_g(dx00,dg2total,dg2,2,order,y,done,msh);
assert(all(done));
end
%%
function [y,done]=combine_g(x0,bd_gval,gval,order,ordlist,y,done,msh)
[isord,loc]=ismember(order,ordlist);
[nw,ngnw,nvec]=deal(length(msh.t),size(gval,1),size(gval,2));
ng=ngnw/nw;
x0=reshape(x0,ng,nvec);
if isord
    done(loc)=true;
    res=msh.wJ*[x0;bd_gval(ng+1:end,:)];
    y([msh.int,msh.id],loc)=...
        {reshape(res,ng,nw,nvec),reshape(gval,ng,nw,nvec)};
end
end