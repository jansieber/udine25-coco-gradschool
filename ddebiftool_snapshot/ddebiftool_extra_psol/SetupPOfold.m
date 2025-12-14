%% SetupPOfold - Initialize continuation of folds of periodic orbits
%%
function varargout=SetupPOfold(funcs,branch,ind,varargin)
%% Inputs
%
%  * |funcs|: functions used for DDE
%  * |branch|: branch of psols along which fold was discovered
%  * |ind|: index in points array that is closest to fold (for initial
%      guess)
%
% optional inputs
%
% * |contpar| (integer default |[]|): set of continuation parameters (if []
%   then free parameters of input branch are used)
% * |sys_deri| (default |1e-4|): used for finite differencing when approximating
%   jacobian of rhs, will be replaced by funcs.sys_deri if funcs.sys_deri is
%   provided by user
% * |correc| (logical, default true): apply |p_correc| to first points on fold
%   branch
% * |dir| (integer, default |[]|): which parameter to vary initially along fold
%   branch (pbranch has only single point if |dir| is empty)
% * |step| (real, default |1e-3|): size of initial step if |dir| is non-empty
% * |hjac| (default |1e-4|) deviation for numerical derivatives if needed
%
% All other named arguments are passed on to pbranch.method.continuation,
% pbranch.method.point and pbranch.parameter
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |pbranch|: fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
varargout=cell(1,3);
[varargout{:}]=SetupPOEV1(funcs,branch,ind,varargin{:});
end
