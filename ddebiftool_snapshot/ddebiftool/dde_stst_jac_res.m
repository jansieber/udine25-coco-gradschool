function varargout=dde_stst_jac_res(funcs,pt,free_par,varargin)
%% residual and jacobian for equilibrium problem
% function [J,res]=dde_stst_jac_res(funcs,pt,free_par,varargin)
% INPUT:
%   funcs problem functions
%	pt current fold solution guess
%	free_par free parameter numbers
% OUTPUT: 
%	J jacobian in R^(n+s x n+p)
%	res residual in R^(n+s x 1)
%
%%
default={'output','J'};
options=dde_set_options(default,varargin(2:end),'pass_on');
[ntaum1,max_rhs_xderiv]=dde_num_delays(funcs);
ntau=ntaum1+1;
[ontau,zerorep]=deal(ones(ntau,1),zeros(size(pt.x,1),ntau*max_rhs_xderiv));
xx=[pt.x(:,ontau),zerorep];
[ind,len]=dde_ind_from_point(pt,free_par);
par=pt.parameter;
res=funcs.wrap_rhs(xx,par);
ieq.de=1:length(res);
varargout=dde_setupoutput('jac_res',options.output,NaN(length(res),len),res,ieq);
if strcmp(options.output,'res')
    return
end
drhs=@(dx,dp)funcs.drhs_mf(xx,par,dx,dp);
J(:,ind.parameter)=drhs(0,      {1,{free_par}});
J(:,ind.x)=    sum(drhs({1,'I'},   0          ), 3);
varargout=dde_setupoutput('jac_res',options.output,J,res,ieq);
end
