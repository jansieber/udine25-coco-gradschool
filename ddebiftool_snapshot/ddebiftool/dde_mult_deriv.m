%% wrapper for arbitrary derivatives of arbitrary order for sys_dirderi, sys_dirdtau
function y=dde_mult_deriv(funcs,type,fcn,directional,varargin)
%% sys_dirdtau: if first extra argument is 'tau' or 'delay'
% cut down x and dx argument format assumptions
[isdelay,iarg]=deal(false,1);
if ismember(type,{'tau','delay'})
    [isdelay,inc0]=deal(true,true);
end
incl_deriv=~isempty(varargin(iarg:end)) && ischar(varargin{iarg}) &&...
    strcmp(varargin{iarg},'incl_deriv');
if incl_deriv
    iarg=iarg+1;
end
if directional
    order=varargin{iarg};
    [x,p]=deal(varargin{iarg+(1:2)});
    if order>0
        dev=reshape(varargin(iarg+(3:4)),2,1);
    else
        dev={};
    end
    opts=varargin(iarg+3+2*double(order>0):end);
else
    %% determine order of derivative and argument format from base argument
    %  order 0=function call
    [x,p,remain]=deal(varargin{iarg+(0:1)},reshape(varargin(iarg+2:end),2,[]));
    optstart=find([cellfun(@ischar,remain(1,:)),true],1,'first');
    order=optstart-1;
    opts=remain(:,optstart:end);
    dev=remain(:,1:order);
end
%% choose function and set dimension of x
xdim_is=2+double(incl_deriv);
if ~isdelay
    nf=size(funcs.lhs_matrixfun(size(x,1)),1);
    xdim_exp=2+double(funcs.delayed_derivs>0);
else % delay
    ndelays=dde_num_delays(funcs);
    nf=ndelays+double(inc0);
    xdim_exp=2;
end
[x,p,dev,xdim]=check_xargs(xdim_exp,xdim_is,x,p,dev);
dims={1,[xdim,2]};
if order>=1 && all(all(cellfun(@(x)isnumeric(x)&&all(x(:)==0),dev),1))
    vecdimx=[size(x,xdim+1:ndims(x)),1];
    vecdimp=[size(p,3:ndims(p)),1];
    nvec=max(length(vecdimx),length(vecdimp));
    vecdimx=[vecdimx,ones(1,max(0,nvec-length(vecdimx)))];
    vecdimp=[vecdimp,ones(1,max(0,nvec-length(vecdimp)))];
    vecdim=max([vecdimx;vecdimp],[],1);
    y=zeros([nf,vecdim]);
    return
end
%% call wrapper for argument expansion
if directional
    y=dir_deriv(fcn{:},dims,order,x,p,dev{:},'nf',nf,opts{:},'splitcomplex',true);
else % and reduction to single direction
    y=mult_deriv(fcn{:},dims,      x,p,dev{:},'nf',nf,opts{:},'splitcomplex',true);
end
end
%%
function [x,p,dev,xdim]=check_xargs(xdim_exp,xdim_is,x,p,dev)
xdim=xdim_exp;
if xdim_exp==xdim_is
    return
end
x=cut_xarg(x,xdim_is);
for i=1:2:numel(dev)
    dev{i}=cut_xarg(dev{i},xdim_is);
end
end
%%
function x=cut_xarg(x,xdim_is)
vecdim=[size(x,xdim_is+1:ndims(x)),1];
nvec=prod(vecdim);
x=reshape(x,size(x,1),size(x,2),size(x,3),nvec);
x=reshape(x(:,:,1,:),[size(x,1),size(x,2),vecdim]);
end