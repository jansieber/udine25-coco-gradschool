%% Initialize continuation of Fold bifurcations
%%
function varargout=SetupFold(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of stst along which Fold was discovered
% * |ind|: index of approximate Fold point
%
% Important optional inputs (name-value pairs)
% outputs
% foldbranch: Fold branch with first point (or two points)
% suc: flag whether corection was successful
%
% Important name-value pair inputs
%
% * |'contpar'| (integers default |[]|): index of additional continuation parameters
%  (in addition to free pars of branch
% * |'correc'| (logical, default true): apply |p_correc| to first points on fold
%   branch
% * |'dir'| (integer, default |[]|): which parameter to vary initially along fold
%   branch (foldbranch has only single point if dir is empty)
% * |'step'| (real, default |1e-3|): size of initial step if dir is non-empty
% * |'outputfuncs'| (default |false|): set to |true| to have 3 outputs,
% with |funcs| first.
%
% All other named arguments are passed on to fields of |foldbranch|.
%% Outputs
% 
% * (if option |'outputfuncs'| is |true|, first output is |funcs|)
% * |foldbranch|: Fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
% Parameter limits for etc are inherited from branch, unless overridden by
% optional input arguments.
%
% $Id: SetupFold.m 309 2018-10-28 19:02:42Z jansieber $
%
%% process options
default={'outputfuncs',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% wrapper around SetupStstBifurcation to ensure backward compatibility
varargout=cell(1,nargout);
[varargout{:}]=SetupStstBifurcation(funcs,branch,ind,'fold',...
    'outputfuncs',options.outputfuncs,pass_on{:});
end
