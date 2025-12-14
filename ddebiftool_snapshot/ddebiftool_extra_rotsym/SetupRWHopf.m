%% SetupRWHopf - Initialize continuation of Hopf bifurcations of relative equilibria
%%
function varargout=SetupRWHopf(funcs,branch,ind,varargin)
%% 
% simple wrapper to have a sensible name (see demo rotsym_demo how to use this function).
default={'outputfuncs',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
out=cell(1,3);
[out{:}]=SetupHopf(funcs,branch,ind,'excludefreqs',0,pass_on{:},'outputfuncs',true);
out=dde_setupoutput('setup',out{:},options.outputfuncs);
varargout=out;
end
