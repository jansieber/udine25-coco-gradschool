%% Initialize continuation of Hopf bifurcations
%%
function varargout=SetupHopf(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of stst along which Hopf was discovered
% * |ind|: index of approximate Hopf point
%
% Important name-value pair inputs
%
% * |'contpar'| (integers default |[]|): index of continuation parameters
%  (replacing free pars of branch)
% * |'correc'| (logical, default true): apply |p_correc| to first points on
% hopf branch
% * |'dir'| (integer, default |[]|): which parameter to vary initially
% along Hopf branch (hbranch has only single point if dir is empty)
% * |'step'| (real, default |1e-3|): size of initial step if dir is non-empty
% * |'outputfuncs'| (default |false|): set to |true| to have 3 outputs,
% with |funcs| first.
%
% All other named arguments are passed on to fields of |hbranch|
%% Outputs
% 
% * (if option |'outputfuncs'| is |true|, first output is |funcs|)
% * |hbranch|: Hopf branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
% Parameter limits for etc are inherited from branch, unless overridden by
% optional input arguments.
%
% $Id: SetupHopf.m 309 2018-10-28 19:02:42Z jansieber $
%
%% wrapper around SetupStstBifurcation to ensure backward compatibility
%branch.method=dde_prepost('add',branch.method,'postprocess','dde_hopf_postprocess',{});
%branch.method=dde_prepost('add',branch.method,'preprocess','dde_jac2square_preprocess',...
%    {'nulldim',1,'free_par',[]});
default={'outputfuncs',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
hopfargs={'point.postprocess','dde_hopf_postprocess'};
%% wrapper around SetupStstBifurcation to ensure backward compatibility
varargout=cell(1,nargout);
[varargout{:}]=SetupStstBifurcation(funcs,branch,ind,'hopf',...
    'outputfuncs',options.outputfuncs,hopfargs{:},pass_on{:});
end
