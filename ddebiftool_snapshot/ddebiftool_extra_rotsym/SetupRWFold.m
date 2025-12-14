%% SetupRWFold - Initialize continuation of folds of relative equilibria
%%
function varargout=SetupRWFold(funcs,branch,ind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of psols along which fold was discovered
% * |ind|: number of point close to fold
%
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |pbranch|: fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%% Optional inputs
%
% * |contpar| (integer default |[]|): set of continuation parameters  
%   (if empty free parameters in argument branch are used)
% * |hbif| (default |1e-3|): used for finite differencing when approximating
%   linearized system,
% * |correc| (logical, default |true|): apply |p_correc| to first points on fold
%   branch
% * |dir| (integer, default |[]|): which parameter to vary initially along fold
%   branch (|pbranch| has only single point if dir is empty)
% * |step| (real, default |1e-3|): size of initial step if dir is non-empty
% * |hjac| (default |1e-6|) deviation for numerical derivatives if needed
%
% all other named arguments are passed on to pbranch.method.continuation,
% pbranch.method.point and pbranch.parameter

%% process options
cond_RWfold=dde_sys_cond_create('name','cond_RWfold',...
    'fun',@(p,pref)sys_cond_RWFold(p,funcs.rotation,pref),...
    'reference',true);
default={'nullparind',branch.parameter.free(end),...
    'usercond',cell(0,1),'initcond',cell(0,1)};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
options.usercond=dde_sys_cond_create(options.usercond);
options.initcond=dde_sys_cond_create(options.initcond);
varargout=cell(1,nargout);
[varargout{:}]=SetupFold(funcs,branch,ind,'outputfuncs',true,...
    'nullparind',options.nullparind,...
    'usercond',[cond_RWfold;options.usercond(:)],...
    'initcond',[cond_RWfold;options.initcond(:)],pass_on{:});
end
