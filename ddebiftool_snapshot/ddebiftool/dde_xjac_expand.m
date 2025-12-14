function [argout,template]=dde_xjac_expand(xpattern,argdims,incl_deriv,maxsize,varargin)
nargs=length(argdims);
argin=reshape(varargin,nargs,[]);
x=argin{1};
other=argin(2:nargs,1);
xdevs=argin(1,2:end);
odevs=argin(2:end,2:end);
nxdims=2+double(incl_deriv);
xdims=size(x,1:nxdims);
xlen=prod(xdims);
[dum,xpatvec]=dde_xpattern('incl_deriv',incl_deriv,'x',x,'xpattern',xpattern); %#ok<ASGLU>
x=reshape(x,[xlen,size(x,nxdims+1:ndims(x)),1]);
is0=cellfun(@(a)isnumeric(a)&&numel(a)==1&&a==0,xdevs);
iscf=cellfun(@(a)iscell(a)&&iscell(a{2}),xdevs);
isid=cellfun(@(a)iscell(a)&&ischar(a{2})&&strcmp(a{2},'I'),xdevs);
isnum=~isid&~iscf&~is0;
xdevs(isnum)=cellfun(@(d){reshape(d,[xlen,size(d,nxdims+1:ndims(d)),1])},xdevs(isnum));
devs=[xdevs;odevs];
zdevs=num2cell(zeros(length(other),1));
lpat=length(xpatvec);
chunksize=ceil(maxsize/numel(x));
argout.chunks=mat2cell(xpatvec,1,diff([0:chunksize:lpat,lpat]));
argout.fun=@(ind)argoutfun(argdims,xdims,x,other,xpatvec,zdevs,devs,nargs,ind);
template=@(fvals)filltemplate(fvals,xlen,xpatvec,xdims);
end
function res=filltemplate(fvals,xlen,xpatvec,xdims)
nf=size(fvals,1);
vecdim=size(fvals,3:ndims(fvals));
nvec=prod(vecdim);
res=zeros(nf,xlen,nvec);
res(:,xpatvec,:)=reshape(fvals,nf,length(xpatvec),nvec);
res=reshape(res,[nf,xdims,vecdim]);
end

function argout = argoutfun(argdims,xdims,x,other,xpatvec,zdevs,devs,nargs,ind)
if nargin<9
    ind=xpatvec;
end
argout=arg_array_expand([1,argdims(2:end)],x,other{:},...
    {1,{ind}},zdevs{:},devs{:});
argout=reshape(argout,nargs,[]);
argout(1,:)=cellfun(@(d){reshape(d,[xdims,size(d,2:ndims(d))])},argout(1,:));
end