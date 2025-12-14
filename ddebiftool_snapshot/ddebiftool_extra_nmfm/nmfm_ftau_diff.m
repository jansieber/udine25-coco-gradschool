%% Differentiate h->sum cj f(xeq+h yj_t) for functional defined by funcs k times
%
% function y=nmfm_ftau_diff(funcs,order,pt,devs,fac)
%
%% Inputs
%
% * funcs: user-provided problem functions, set with set_funcs
% * order: order of derivative to be taken
% * pt: point of kind 'stst' (or local bifurcation), pt.x and pt.parameter
% are needed
% * devs: 1 x ndev array of history functions yj
% * fac; 1 x ndev array of complex numbers cj
%
% devs and fac are obtained from polarization identity in nmfm_mfderiv
% (nmfm_dev_group, nmfm_pol_from_devs) to compute an arbitrary mixed
% derivative.
%
%% Outputs
%
% * y: directional derivative d/dh [sum cj f(x0+h yj)] at h=0
% * y2: lower order estimate (when using finite differences)
% * h: interpolation point spacing used for finte difference approximation
%
% $Id: nmfm_ftau_diff.m 309 2018-10-28 19:02:42Z jansieber $
%%
function [y,y2,h]=nmfm_ftau_diff(funcs,order,pt,devs,fac,varargin)
default={'debug_derivs',0};
options=dde_set_options(default,varargin,'pass_on');
ndev=length(fac);
ndim=size(pt.x,1);
p0=pt.parameter(1,:,ones(1,ndev));
max_rhs_xderiv=dde_num_delays(funcs,'max_rhs_xderiv');
taus=dde_stst_delays(funcs,pt);
ntaus=length(taus);
apply=@(f,x)cell2mat(arrayfun(f,x,'UniformOutput',false));
%% compute derivatives of deviations as needed
max_order=double(logical(funcs.tp_del))*order+max_rhs_xderiv+1;
devs=repmat(devs(:),1,max_order);
for i=1:max_order-1
    devs(:,i+1)=arrayfun(@nmfm_dev_deriv,devs(:,i));
end
xx0vec=repmat(pt.x,1,ntaus,ndev,max_order);
xx0vec(:,:,:,2:end)=0;
devvals=apply(@(d)real(nmfm_dev_call(d,-taus)),reshape(devs(:),1,1,[]));
devvals=reshape(devvals,[],ntaus,ndev,max_order);
if ~funcs.tp_del
    devsel=devvals(:,:,:,1:max_order);
    xdevs=devsel(1:ndim,:,:);
    pdevs=reshape(devvals(ndim+1:end,1,:),1,[],ndev);
    xx0sel=xx0vec(:,:,:,1:max_order);
    yi=funcs.drhs_dir(order,xx0sel,p0,xdevs,pdevs);
else
    devsel=apply(@(i)devvals(:,:,:,1+i+(0:max_rhs_xderiv)),reshape(0:order,1,1,1,[]));
    xdevs=devsel(1:ndim,:,:,:);%nx,ntaup1,ndev,deriv
    pdevs=reshape(devvals(ndim+1:end,1,:,1),[],ndev);
    xx0sel=apply(@(i)xx0vec(:,:,:,1+i+(0:max_rhs_xderiv)),reshape(0:order,1,1,1,[]));    
    yi=nmfm_dirmf_combined(funcs,xx0sel,p0,order,taus,xdevs,pdevs,varargin{:});
    if options.debug_derivs>0
        for i=ndev:-1:1
            f0=@(h)nmfm_ftau_dev(funcs,pt,devs(i),1,h);
            yn(:,i)=num_scalarderiv(f0,order);
        end
        assert(all(abs(yi(:)-yn(:))<eps^(1./order/2)));
    end
end
facvec=repmat(reshape(fac,1,ndev),size(yi,1),1);
y=sum(yi.*facvec,2);
y2=y;
h=0;
end