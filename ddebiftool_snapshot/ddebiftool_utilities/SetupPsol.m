%% Start branch of periodic orbits from Hopf point or bifurcation of periodic orbit branch
%%
function varargout=SetupPsol(funcs,branch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user
% * |branch|: branch of |'stst'| steady state solutions, |'hopf'|
% solutions from which one wants to branch off, or |'psol'| solutions that are bifurcations of psols
% * |ind|: index in |'point'| field of |branch| where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% (for branching off from Hopf)
%
% * |'degree'|: degree used for collocation
% * |'intervals'|: number of collocation intervals (overall mesh size is
% degree x intervals + 1
% * |'hopfcorrection'|: whether Hopf point still needs to be corrected
%
% (general)
%
% * |'contpar'|: index of continuation parameters if different from
% |'branch'|, (this should be set manually when the original branch was a
% two-parameter branch)
% * |'corpar'|: parameters left free for initial correction (if different
% from |contpars|)
% * |'outputfuncs'| (default |false|): set to |true| to have 3 outputs,
% with |funcs| first.
% * |'outputfuncs'| (default |false|): set to |true| to have 3 outputs,
% with |funcs| first.
%
% All other optional inputs are passed on to fields of out branch |per|
%% Outputs
%
% * (if option |'outputfuncs'| is |true|, first output is |funcs|)
% * |per|: branch of periodic orbits with desired settings and two initial
% corrected points
% * |suc|: flag indicating success
% 
% $Id: SetupPsol.m 309 2018-10-28 19:02:42Z jansieber $
%
%%
default={'contpar',[],'corpar',[],'outputfuncs',false,'branch_off',''};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if isempty(options.contpar)
    options.contpar=branch.parameter.free;
end
if isempty(options.corpar)
    % assume that initial correction have to be made in all continuation
    % parameters (to accomodate step condition)
    options.corpar=options.contpar;
end
if isempty(options.branch_off)
    options.branch_off=branch.point(ind).kind;
end
[pfuncs,per,suc]=dde_apply({'SetupPsolFrom_',options.branch_off,''},funcs,branch,ind,...
    'contpar',options.contpar,'corpar',options.corpar,...
    'outputfuncs',true,pass_on{:});
varargout=dde_setupoutput('setup',pfuncs,per,suc,options.outputfuncs);
end
