%% wrapper for directional derivatives of arbitrary order
function y=dfdir_wrap(fun,fmt,u0,udev,dims,argdim)
[n,nvec]=size(u0);
udev=reshape(udev,n*fmt.deg,nvec);
if fmt.deg==0
    y=fmt.fuwrap(fun,0,argdim,u0);
    return
end
cf=(dec2bin(2^(fmt.deg-1):2^fmt.deg-1)-'0')*2-1;
nfac=size(cf,1);
lfac=prod(cf,2)/factorial(fmt.deg)/2^(fmt.deg-1);
dev=kron(cf,eye(n))*udev;
dev=reshape(dev,n,nfac*nvec);
urep=reshape(repmat(u0,nfac,1),n,nfac*nvec);
yfac=fmt.fuwrap(fun,fmt.deg,argdim,urep,dev);
nf=size(yfac,1);
yfac=reshape(yfac,nf*nfac,nvec);
y=kron(lfac',eye(nf))*yfac;
if length(dims)>1
    y=full(y);
end
y=reshape(y,[nf,dims]);
end
