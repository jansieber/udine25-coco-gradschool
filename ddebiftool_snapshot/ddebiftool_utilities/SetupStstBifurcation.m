%% Initialize continuation of equilibrium bifurcation (fold or Hopf)
function varargout=SetupStstBifurcation(funcs,branch,ind,kind,varargin)
%% Inputs
%
% * |funcs|: functions used for DDE
% * |branch|: branch of stst along which bifurcation was discovered
% * |ind|: index of approximate bifurcation point
% * |type|: tpye of bifurcation, string: at the moment either |'fold'| or |'hopf'|
%
% Important name-value pair inputs
%
% * |'contpar'| (integers default |[]|): index of continuation parameters
%  (replacing free pars of branch)
% * |'correc'| (logical, default true): apply |p_correc| to first points on
% bifurcation branch
% * |'dir'| (integer, default |[]|): which parameter to vary initially
% along bifurcation branch (bifbranch has only single point if dir is empty)
% * |'step'| (real, default |1e-3|): size of initial step if dir is non-empty
% * |'outputfuncs'| (default |false|): set to |true| to have 3 outputs,
% with |funcs| first.
%
% All other named arguments are passed on to fields of |bifbranch|
%% Outputs
% 
% * |bifbranch|: branch of bifurcation points with first point (or two points)
% * |suc|: flag whether corection was successful
% * |funcs|: same as |funcs|, unless branch point
% * (if option |'outputfuncs'| is |true|, first output is |funcs|)
%
% Parameter limits for bifbranch etc are inherited from branch, unless overridden by
% optional input arguments.
%
% $Id: SetupStstBifurcation.m 309 2018-10-28 19:02:42Z jansieber $
%
%% process options
default={'contpar',[],'correc',true,'dir',[],'step',1e-3,...
    'usercond',cell(0,1),'initcond',cell(0,1),...
    'outputfuncs',false,'nullparind',zeros(0,1),'output','point'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
point=branch.point(ind);
if ~isfield(point,'stability') || isempty(point.stability)
    point.stability=p_stabil(funcs,point,branch.method.stability);
end
%% check if extra parameters are needed
nnull=length(options.nullparind);
if nnull>0
    %% extend problem
    [funcs,extra_freepar,initfuncs]=feval(['set_',kind,'funcs'],...
        funcs,point,'nullparind',options.nullparind,...
        'usercond',options.usercond,'initcond',options.initcond);
    options.contpar=[options.contpar,extra_freepar];
    point.parameter(funcs.ip.nullparind(:,2))=0;
else
    %% add new sys_cond if given by user, usercond instead of initcond
    initfuncs=dde_add_cond('SetupStstBifurcation',funcs,options,point,'initcond');
    funcs=dde_add_cond('SetupStstBifurcation',funcs,options,point,'usercond');
end
%% initialize branch of bifurcations (bifbranch)
bifbranch=branch;
bifbranch=replace_branch_pars(bifbranch,options.contpar,[{'branchtype',kind},pass_on]);
%% find critical frequency (for Hopf)
[freq,point]=dde_stst_critical_freq(initfuncs,point,kind,branch.method.stability,pass_on{:});
%% create initial guess for correction
[pini0,v,sv]=dde_stst_critical_space(initfuncs,point,kind,freq,bifbranch.method.point,...
    pass_on{:},'output','point');
if options.outputfuncs
    pini0.nvec=[];
    pini0.nmfm=[];
end
bifbranch.point=pini0;
if strcmp(options.output,'sv')
    varargout={sv(end:-1:1),v,funcs};
    return
end
%% correct and add 2nd point if desired
[bifbranch,suc]=correct_ini(funcs,bifbranch,pini0,...
    options.dir,options.step,options.correc);
varargout=dde_setupoutput('setup',funcs,bifbranch,suc,options.outputfuncs);
end
