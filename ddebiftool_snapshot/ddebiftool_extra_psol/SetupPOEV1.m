%% SetupPOEV1 - Initialize continuation of periodic orbits w exrta FM=1
%%
function varargout=SetupPOEV1(funcs,branch,ind,varargin)
%% Inputs
%
%  * |funcs|: functions used for DDE
%  * |branch|: branch of psols along which EV1 was discovered
%  * |ind|: index in points array that is closest to EV1 (for initial
%      guess)
%
% optional inputs
%
% * |contpar| (integer default |[]|): set of continuation parameters (if []
%   then free parameters of input branch are used)
% * |sys_deri| (default |1e-4|): used for finite differencing when approximating
%   jacobian of rhs, will be replaced by funcs.sys_deri if funcs.sys_deri is
%   provided by user
% * |correc| (logical, default true): apply |p_correc| to first points on
%    POEV1 branch
% * |dir| (integer, default |[]|): which parameter to vary initially along
%    POEV1 branch (pbranch has only single point if |dir| is empty)
% * |step| (real, default |1e-3|): size of initial step if |dir| is non-empty
% * |hjac| (default |1e-4|) deviation for numerical derivatives if needed
% * |'POEV1_norm'| (default true) add extra condition that eigenvectors
% has length 1
% * |'POEV1_phase_condition'| (default true) add extra condition that eigenvectors
% has length 1
%
% All other named arguments are passed on to pbranch.method.continuation,
% pbranch.method.point and pbranch.parameter
%% Outputs
%
% * |pfuncs|: functions used for extended DDE
% * |pbranch|: POEV1 branch with first point (or two points)
% * |suc|: flag whether corection was successful
%
%% process options
default={'contpar',[],'correc',true,'dir',[],'step',1e-3,'output','funcs'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
branch.point=branch.point(ind); % remove all points but approx solution with FM=1
%% If branch is psol bifurcation remove extra args
[funcs,branch]=PsolFromPsolbif(funcs,branch);
%% initialize branch of POEV1 solutions (pbranch) and pass on optional args
pbranch=replace_branch_pars(branch,options.contpar,pass_on);
point=dde_psol_create(pbranch.point);
[pfuncs,extra_freepar,initfuncs]=set_POfoldfuncs(funcs,point,...
    branch.method,pass_on{:});
%% required amendments of structures for extended system
free_par=[pbranch.parameter.free,extra_freepar];
[dum,is]=unique(free_par); %#ok<ASGLU>
pbranch.parameter.free=free_par(sort(is));
pbranch.method.point.extra_condition=1;
pbranch.method.point.phase_condition=0;
%% create initial guess for correction
out=cell(1,3);
[out{:}]=POEV1Init(initfuncs,point,pbranch.method,'output',options.output,pass_on{:});
if strcmp(options.output,'sv')
    varargout=[out(1:2),{pfuncs}];
    return
end
pfoldini0=out{1};
%% correct initial guess and find 2nd point along branch if desired
[pbranch,suc]=correct_ini(pfuncs,pbranch,pfoldini0,...
    options.dir,options.step,options.correc);
varargout={pfuncs,pbranch,suc};
end
%
