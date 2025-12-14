function [A,taus,max_rhs_xderiv]=dde_stst_linearize_rhs(funcs,point,varargin)
%% linearize r.h.s. in stst
%
%%
default={'deriv',{{1,'I'},0},'repeat',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
max_rhs_xderiv=dde_num_delays(funcs,'max_rhs_xderiv');
[taus,xx]=dde_stst_delays(funcs,point,'repeat',true,'max_xderiv',max_rhs_xderiv,pass_on{:});
A = funcs.drhs_mf(xx,point.parameter,options.deriv{:});
end
%%
