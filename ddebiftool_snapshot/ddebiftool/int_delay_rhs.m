function [res,gvals,sm,msh]=int_delay_rhs(gfun,msh,order,xx,p,bd,val,dxx,dp,dbd,dval)
%% approximate algebraic equation for integral distributed delay
%
% 0=sum wt(i)*tau*g(t(i)*tau,xx(:,t-t(i)*tau),p)-val
%
% where wt(i) and t(i) are predefined integration weights and nodes on the
% interval [0,1], if the interval is finite, or on [0,infty] if the
% interval is infinite (Laguerre or generalised Laguerre approximation)
%%
[ng,nw,nx,nvec]=deal(size(val,1),length(msh.w),size(xx,1),size(xx,3));
bd=reshape(bd,1,nvec);
wt=repmat(reshape(msh.w(:)*bd,1,nw,nvec),[ng,1,1]);
s=reshape(msh.t(:)*bd,1,1,nw*nvec);
if nargout>=3
    sm=reshape(s,nw,nvec);
end
xs=reshape(xx,nx,1,nw*nvec);
ps=reshape(repmat(p,1,nw,1),1,size(p,2),nw*nvec);
gvals=reshape(gfun(0,s,xs,ps),ng,nw,nvec);
if order==0
    res=reshape(-val(:,1,:)+sum(wt.*gvals,2),ng,nvec);
    return
end
dbd=reshape(dbd,1,nvec);
dwt=reshape(repmat(reshape(msh.w(:)*dbd,1,nw*nvec),[ng,1,1]),ng,nw,nvec);
ds=reshape(msh.t(:)*dbd,1,1,nw*nvec);
dxs=reshape(dxx,nx,1,nw*nvec);
dps=reshape(repmat(dp,1,nw,1),1,size(dp,2),nw*nvec);
dgvals=reshape(gfun(order,s,xs,ps,ds,dxs,dps),ng,nw,nvec);
if order>1
    dgm1vals=reshape(gfun(order-1,s,xs,ps,ds,dxs,dps),ng,nw,nvec);
else
    dgm1vals=gvals;
end
dgtotal=order*dwt.*dgm1vals+wt.*dgvals;
res=reshape(sum(dgtotal,2),ng,nvec);
if order==1
    res=res-reshape(dval(:,1,:),ng,nvec);
end
end
