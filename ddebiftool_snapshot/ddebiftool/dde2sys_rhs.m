function y=dde2sys_rhs(deg,ix,ip,f,fargs,xx,p,dx,dp)
%% reshape argument format of f to the one used by sys_rhs and sys_dirderi
[nx,nd]=size(xx,1:2);
np=size(p,2);
vecdim=[size(xx,3:ndims(xx)),1];
nvec=prod(vecdim);
xx=reshape(xx,nx,nd,nvec);
xc=xc_from_xx(xx,ix);
p=reshape(p,np,nvec);
p=p(ip,:);
pc=num2cell(p,2);
if deg==0
    y=f(fargs{:},xc{:},pc{:});
    ny=size(y,1);
    y=reshape(y,[ny,vecdim]);
    return
end
dx=reshape(dx,nx,nd,nvec);
dxc=xc_from_xx(dx,ix);
dp=reshape(dp,np,nvec);
dp=dp(ip,:);
dpc=num2cell(dp,2);
y=f(fargs{:},xc{:},pc{:},dxc{:},dpc{:});
ny=size(y,1);
y=reshape(y,[ny,vecdim]);
end
%%
function xc=xc_from_xx(xx,ix)
if isnumeric(ix)
    xx=xx(ix,:,:);
    xc=num2cell(permute(xx,[2,3,1]),[1,2]);
else
    xc=num2cell(permute(xx,[1,3,2]),[1,2]);
end
end