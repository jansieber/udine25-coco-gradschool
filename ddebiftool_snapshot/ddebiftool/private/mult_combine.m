function y=mult_combine(fundirderiv,order,dims,matmul,x,dx,groups)
if isnumeric(dims)
    dims={1,dims};
end
[outdims,argdims]=deal(dims{:});
[coeffs,factors,used]=polarization_coeffs(order,groups);
vecdims=[size(x{1},argdims(1)+1:ndims(x{1})),1];
fmt=cellfun(@(a,f){size(a,1:f)},x(:),num2cell(argdims(:)));
nvec=prod(vecdims);
nfac=length(factors);
ndxcols=size(dx,2);
alen=cellfun(@prod,fmt);
devvec=cellfun(@(d){reshape(d,[],nvec)},dx);
devmat=reshape(permute(reshape(cat(1,devvec{:}),[],ndxcols,nvec),[2,1,3]),ndxcols,[]); % order x sum(prod(fmt{:}))*nvec
dev=reshape(permute(reshape(matmul(coeffs,devmat(used,:)),nfac,[],nvec),[2,1,3]),[],nfac*nvec);
if isempty(dev)
    y=fundirderiv(x{:},dx{:,1:min(1,ndxcols)});
    return
end
devc=mat2cell(dev,alen,nfac*nvec);
devfmt=cellfun(@(a,f){reshape(a,[f,nfac*nvec,1])},devc,fmt);
baserep=cellfun(@(a,f){reshape(repmat(reshape(a,[],1,nvec),1,nfac,1),...
    [f,nfac*nvec,1])},x,fmt);
yfac=fundirderiv(baserep{:},devfmt{:});
nf=size(yfac,1:outdims);
numf=prod(nf);
yfac=reshape(yfac,numf*nfac,nvec);
y=kron(factors.',eye(numf))*yfac;
if length(vecdims)>1
    y=full(y);
end
y=reshape(y,[nf,vecdims]);
end