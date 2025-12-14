function [taus,xtau_values,pvec]=nmfm_dev_delays(funcs,devs,nx,varargin)
%% find delays for equilibria plus small deviations: 
% devs is an array of deviation structures (numel nvec), including
% parameters 
% 
% output are (ntau+1) x nvec array taus giving the value of tau(k,dev(j))
% max_rhs_xderiv as stored in funcs, xtau_values is 
% n x (ntau+1) x nvec x (maxorder+1) array x(i,j,k,l) containing
% x_i^(l)(-tau_j) for dev(k)
ntaum1=dde_num_delays(funcs); 
default={'d_nr',0:ntaum1,'max_xderiv',1};
options=dde_set_options(default,varargin,'pass_on');
max_order=options.max_xderiv;
np=size(devs(1).v,1)-nx;
nvec=numel(devs);
devs=repmat(devs(:),1,max_order+1);
%% differentiate all deviations up to maxorder
for i=1:max_order
    devs(:,i+1)=arrayfun(@nmfm_dev_deriv,devs(:,i));
end
apply=@(f,a,rs)cell2mat(arrayfun(f,reshape(a,rs{:}),'UniformOutput',false));
pvec=apply(@(d)nmfm_dev_call(d,0,0,nx+(1:np)),devs(:,1),{1,[]});
pvec=reshape(pvec,1,np,nvec);
xhistfun=@(it,t,deriv)eval_devs(it,t,deriv,devs,nx);
[taus,xtau_values]=dde_fun_delays(funcs,xhistfun,...
    nvec,pvec,'d_nr',options.d_nr,...
    'max_xderiv',options.max_xderiv);
end
%%
function x=eval_devs(itau,negtau,derivs,devs,nx)
apply=@(f,a,rs)cell2mat(arrayfun(f,reshape(a,rs{:}),'UniformOutput',false));
ntauc=length(itau);
negtau=reshape(negtau,ntauc,[]);
nvec=size(negtau,2);
for k=length(derivs):-1:1
    x(:,:,:,k)=apply(@(i)nmfm_dev_call(devs(i,derivs(k)+1),negtau(:,i),0,1:nx),1:nvec,{1,1,[]});
end
end
