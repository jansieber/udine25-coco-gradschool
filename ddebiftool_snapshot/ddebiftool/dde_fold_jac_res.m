function varargout=dde_fold_jac_res(funcs,pt,free_par,method,varargin)
%% Jacobian and residual of nonlinear system for fold
% function [J,res,ieq]=dde_fold_jac_res(funcs,pt,free_par,...)
% INPUT:
%   funcs problem functions
%	pt current fold solution guess
%	free_par free parameter numbers
% OUTPUT: 
%	J jacobian in R^(2n+1+s x 2n+p)
%	res residual in R^(2n+1+s x 1)
%   ieq index vectors for equations 
%
%%
default={'output','J'};
options=dde_set_options(default,varargin,'pass_on');
[ind,len]=dde_ind_from_point(pt,free_par);
[ntaum1,max_rhs_xderiv]=dde_num_delays(funcs);
ntau=ntaum1+1;
[ontau,zerorep]=deal(ones(ntau,1),zeros(size(pt.x,1),ntau*max_rhs_xderiv));
xx=[pt.x(:,ontau),zerorep];
vx=[pt.v(:,ontau),zerorep];
[ifree,ifree_ext,nuserpar]=get_ifree(funcs,pt,free_par);
par=pt.parameter(1:nuserpar);
vpar=0*par;
vpar(ifree)=pt.parameter(ifree_ext);
%% residual
r0=funcs.wrap_rhs(xx,par);
r1=funcs.drhs_dir(1,xx,par,vx,vpar);
res=[r0;r1];
ieq.de=1:length(r0);
ieq.dvar=ieq.de(end)+ieq.de;
%% append normalization condition
[resnorm,Jcnorm,ieqnorm]=dde_stst_nullnorm_cond(pt,ifree_ext,method);
res=[res;resnorm];
ieq.norm=ieq.dvar(end)+ieqnorm;
varargout=dde_setupoutput('jac_res',options.output,NaN(length(res),len),res,ieq);
if strcmp(options.output,'res')
    return
end
%% Jacobian
[i_noext,ip_noext]=setdiff(free_par,ifree_ext);
[dum,ip_ext]=ismember(ifree_ext,free_par); %#ok<ASGLU>
drhs=@(varargin)funcs.drhs_mf(xx,par,varargin{:});
J=zeros(2*length(ieq.de),len);
J(ieq.de,ind.parameter(ip_noext(:)))=drhs(0,{1,{i_noext}});
J(ieq.de,ind.x)=sum(drhs({1,'I'},0),3);
J(ieq.dvar,ind.x)=sum(drhs(vx,vpar,{1,'I'},0),3);
J(ieq.dvar,ind.v)=J(ieq.de,ind.x);
J(ieq.dvar,ind.parameter([ip_noext(:);ip_ext(:)]))=...
    [drhs(vx,vpar,0,{1,{i_noext}}),drhs(0,{1,{ifree}})];
%% append normalization condition
J_cond=dde_x_from_point(Jcnorm(:),free_par);
J=[J;J_cond.'];
varargout=dde_setupoutput('jac_res',options.output,J,res,ieq);
end
%% find out which free parameters were part of the basic problem
% ifree is part of basic zero problem (eg, the rotation speed for
% rotational symmetry), ifree_ext is its linear deviation, iext are other
% parameters (for continuation)
function  [ifree,ifree_ext,nuserpar]=get_ifree(funcs,pt,free_par)
[ifree,ifree_ext,nuserpar]=deal(zeros(1,0),zeros(1,0),length(pt.parameter));
if ~isfield(funcs,'ip') || ~isfield(funcs.ip,'nullparind')
    return
end
ifree=intersect(funcs.ip.nullparind(:,1)',free_par);
ifree_ext=intersect(reshape(funcs.ip.nullparind(:,2:end),1,[]),free_par);
nuserpar=funcs.ip.nuserpar;
end
