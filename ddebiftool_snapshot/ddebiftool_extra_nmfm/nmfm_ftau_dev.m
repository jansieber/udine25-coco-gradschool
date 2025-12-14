%% Evaluate function h->sum cj f(xeq+h yj_t) for functional defined by funcs
%
% function y=nmfm_ftau_dev(funcs,pt,devs,fac,h)
%
%% Inputs 
%
% * funcs: structure containing problem functions
% * pt: constant point to which deviations are added
% * devs (1 x ndev): array of deviation structures
% * fac (1 x ndev): factors cj in weighted sum
% * h (1 x nh): distance of deviation (vectorized)
%
%% Outputs
%
% * y: (ny x nh) array sum cj f(pt.x+h dev(j),pt.parameter)
%
% $Id: nmfm_ftau_dev.m 309 2018-10-28 19:02:42Z jansieber $
%%
function y=nmfm_ftau_dev(funcs,pt,devs,fac,h)
ndev=length(fac);
nh=length(h);
ndim=length(pt.x);
[ntau,max_order]=dde_num_delays(funcs); %#ok<ASGLU>
%% generate deviations from equilibrium
eqfun=nmfm_dev_fun([pt.x;pt.parameter(:)]);
[hm,dm]=ndgrid(h,devs(:));
devs_h=arrayfun(@(harg,d)nmfm_dev_ax([1,harg],[eqfun,d]),hm,dm);
%% differentiate deviations
% after this devs_h(i,j,k) is x0+h(i)*diff(devs(j),k-1)
devs_h=repmat(devs_h,1,1,max_order+1);
for i=1:max_order
    devs_h(:,:,i+1)=arrayfun(@nmfm_dev_deriv,devs_h(:,:,i));
end
%% find values of history where to compute functional
% done differently if delays are state-dependent or constant
[taus,xtau_values,pvec]=nmfm_dev_delays(funcs,devs_h(:,:,1),ndim,'max_xderiv',max_order); %#ok<ASGLU>
%% evaluate functional sys_rhs in all deviations
yj=funcs.wrap_rhs(xtau_values,pvec);%NaN(ndim,nvec);
%% add up functional values, weighted with factors fac (cj)
yj=reshape(yj,[],nh,ndev);
fac=reshape(fac,[1,1,ndev]);
facvec=repmat(fac,[size(yj,1),nh,1]);
y=sum(yj.*facvec,3);
end
