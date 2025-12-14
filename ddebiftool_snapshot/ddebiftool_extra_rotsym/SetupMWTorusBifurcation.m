%% SetupMWTorusBifurcation - Initialize continuation of torus bifurcations of relative periodic orbits
%%
function [trfuncs,trbranch,suc]=SetupMWTorusBifurcation(funcs,branch,ind,varargin)
%%
% simple wrapper to have a sensible name and set number of trivial Floquet
% multipliers to two (see demo rotsym_demo how to do this).
%
% $Id: SetupMWTorusBifurcation.m 319 2019-01-31 02:14:56Z jansieber $
%
default={'nullparind',branch.parameter.free(end),'usercond',cell(0,1)};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[trfuncs,trbranch,suc]=SetupTorusBifurcation(funcs,branch,ind,'nremove',[1,1],...
    'nullparind',options.nullparind,...
    'usercond',options.usercond,pass_on{:});
end
