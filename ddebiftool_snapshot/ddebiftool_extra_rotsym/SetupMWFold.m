%% SetupMWFold - Initialize continuation of folds of modulated waves
%%
function varargout=SetupMWFold(funcs,branch,ind,varargin)
%% Inputs
% 
% * |funcs|: functions used for rotationally symmetric DDE
% * |branch|: branch of psols along which fold was discovered
% * |ind| number of point close to fold
%
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |pbranch|: fold branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%% Optional inputs
% 
% * |contpar| (integers default |[]|): index of continuation parameters 
%   (replacing free parameters in argument branch)
% * |hbif| (default |1e-3|): used for finite differencing when approximating
%   linearized system, replace by |funcs.sys_deri| if |funcs.sys_deri| is
%   analytical
% * |correc| (logical, default |true|): apply |p_correc| to first points on fold
%   branch
% * |dir| (integer, default |[]|): which parameter to vary initially along fold
%   branch (|pbranch| has only single point if |dir| is empty)
% * |step| (real, default |1e-3|): size of initial step if dir is non-empty
%
% all other named arguments are passed on to pbranch.method.continuation,
% pbranch.method.point and pbranch.parameter
%
% $Id: SetupMWFold.m 374 2019-09-14 14:02:58Z jansieber $
%

%% process options
branch.point=branch.point(ind);
[funcs,branch]=PsolFromPsolbif(funcs,branch);
cond_MWfold=dde_sys_cond_create('name','cond_MWFold','fun',...
    @(p,pref)sys_cond_MWFold(p,pref,funcs),'reference',true);
default={'nullparind',branch.parameter.free(end),...
    'usercond',cell(0,1),'initcond',cell(0,1)};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
options.usercond=dde_sys_cond_create(options.usercond);
options.initcond=dde_sys_cond_create(options.initcond);
varargout=cell(1,nargout);
[varargout{:}]=SetupPOfold(funcs,branch,1,...
    'nullparind',options.nullparind,...
    'usercond',[cond_MWfold;options.usercond(:)],...
    'initcond',[cond_MWfold;options.initcond(:)],pass_on{:});
end
