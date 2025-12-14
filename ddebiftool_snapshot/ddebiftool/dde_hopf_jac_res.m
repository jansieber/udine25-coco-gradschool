function varargout=dde_hopf_jac_res(funcs,pt,free_par,method,varargin)
%% Jacobian for Hopf problem
% function [J,res]=dde_hopf_jac_res(funcs,pt,free_par,method,pref,varargin)
% INPUT:
%   funcs problem function
%	x current Hopf solution guess in R^n
%	omega current Hopf frequency guess in R
%	v current eigenvector guess in C^n
%	par current parameters values in R^p
%	free_par free parameter numbers in N^d 
%	c normalization vector of v in C^(1 x n)
% OUTPUT: 
%	J jacobian in R^(3n+2+s x 3n+1+p)
%	res residual in R^(3n+2+s x 1)
%   ieq index vectors for equations 
%
%%
default={'pref',repmat(pt,0,1),'output','J'};
options=dde_set_options(default,varargin,'pass_on');
pref=options.pref;
[ind,len]=dde_ind_from_point(pt,free_par);
%% residual
[Jstst,rstst,ieq]=dde_stst_jac_res(funcs,pt,free_par);
[lambda,v]=deal(1i*pt.omega,pt.v);
Delta=@(varargin)ch_matrix(funcs,pt.x,pt.parameter,lambda,varargin{:});
Delta0=Delta();
res=[...
    rstst;...
    real(Delta0*v);
    imag(Delta0*v)];
[ieq.dvar.re,ieq.dvar.im]=deal(ieq.de(end)+ieq.de,2*ieq.de(end)+ieq.de);
%% append normalization condition if required by method
[resnorm,Jcnorm,ieqnorm]=dde_stst_nullnorm_cond(pt,[],method);
ieq.norm=ieq.dvar.im(end)+ieqnorm;
ieqlen=ieq.dvar.im(end)+length(ieqnorm);
%% if pref is non-empty minimize phase difference
[resph,Jph,ieqph]=dde_stst_nullphase_cond(pt,pref,[],method);
ieq.phase=ieqlen+ieqph;
%% append extra conditions
res=[res;resnorm;resph];
varargout=dde_setupoutput('jac_res',options.output,NaN(length(res),len),res,ieq);
if strcmp(options.output,'res')
    return
end
%% Jacobian
dDdxv=Delta('dx',v,'devx',{1,'I'},'order',1);
dDdpv=Delta('dx',v,'devp',{1,{free_par}},'order',1);
dDvom=Delta('dx',v,'deri',1)*1i;
%% assemble jacobian
J=zeros(ieq.dvar.im(end),len);
J(ieq.de,[ind.x(:);ind.parameter(:)])=Jstst(:,[ind.x(:);ind.parameter(:)]);
ic=                 [ind.x(:);   ind.v.re(:); ind.v.im(:); ind.omega;  ind.parameter(:)];
J(ieq.dvar.re,  ic)=[real(dDdxv),real(Delta0),-imag(Delta0), real(dDvom),real(dDdpv)];
J(ieq.dvar.im,  ic)=[imag(dDdxv),imag(Delta0), real(Delta0), imag(dDvom),imag(dDdpv)];
J_cond=dde_x_from_point([Jcnorm(:);Jph],free_par);
J=[J;J_cond.'];
varargout=dde_setupoutput('jac_res',options.output,J,res,ieq);
end
