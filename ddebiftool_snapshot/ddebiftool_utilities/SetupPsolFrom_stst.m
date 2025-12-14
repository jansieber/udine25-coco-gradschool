%% Branch off at Hopf point (point is either of kind hopf or stst)
%%
function varargout=SetupPsolFrom_stst(funcs,ststbranch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user
% * |ststbranch|: branch of |'stst'| steady state solutions or |'hopf'|
% solutions from which one wants to branch off
% * |ind|: index in |'point'| field of |ststbranch| which is close to Hopf
% point where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |'degree'|: degree used for collocation
% * |'intervals'|: number of collocation intervals (overall mesh size is
% degree x intervals + 1
% * |'hopfcorrection'|: whether Hopf point still needs to be corrected
% * |'contpar'|: index of continuation parameters if different from
% |'ststbranch'|
% * |'corpar'|: parameters left free for initial correction (if different
% from |contpars|)
%
% All other optional inputs are passed on to fields of out branch |per|
%% Outputs
%
% * (|pfuncs|: possibly modified rhs functions. This output is given only
% if requested by optional argument |'outputfuncs', true|.)
% * |per|: branch of periodic orbits with desired settings and two initial
% corrected points
% * |suc|: flag indicating success
% 
% $Id: SetupPsolFrom_stst.m 369 2019-08-27 00:07:02Z jansieber $
%
%%
default={'radius',0.01,'contpar',[],'corpar',[],...
    'includehopf',true,...
    'degree',3,'intervals',20,'usercond',[],'outputfuncs',false};
[hopfargs,args]=dde_options_filter('SetupHopf',varargin);
[options,pass_on]=dde_set_options(default,args,'pass_on');
% create branch per of periodic solutions starting from an
% approximate Hopf point num on a branch br of steady states
if isempty(options.contpar)
    options.contpar=ststbranch.parameter.free;
end
if isempty(options.corpar)
    % assume that initial correction have to be made in all continuation
    % parameters (to accomodate step condition)
    options.corpar=options.contpar;
end
[hfuncs,hbr,suc]=SetupHopf(funcs,ststbranch,ind,'outputfuncs',true,...
    'contpar',options.contpar,'corpar',options.corpar,...
    'includehopf',options.includehopf,pass_on{:},hopfargs{:});
if suc==0
    varargout=dde_setupoutput('setup',hfuncs,[],suc,options.outputfuncs);
    return
end
%% set method and branch parameters
per=df_brnch(funcs,options.contpar,'psol');
per.parameter=ststbranch.parameter;
per.parameter.free=options.contpar;
per=replace_branch_pars(per,options.contpar,pass_on);
mth=per.method.point;
mthargs={'collocation_parameters',mth.collocation_parameters};
%% generate first initial solutions
hopf=hbr.point(1);
deg_psol=p_topsol(funcs,hopf,0,options.degree,options.intervals,...
    pass_on{:},mthargs{:});
[dev_psol,step_cond]=p_topsol(funcs,hopf,...
    options.radius,options.degree,options.intervals,pass_on{:});
funcs=dde_add_cond('SetupPsolFrom_stst',funcs,options,dev_psol);
[psol,suc]=p_correc(funcs,dev_psol,options.contpar,step_cond,per.method.point,0, ...
    dev_psol,pass_on{:});
if suc==0
    varargout=dde_setupoutput('setup',funcs,[],suc,options.outputfuncs);
    return;
end
per.point=deg_psol;
per.point.profile=repmat(hopf.x,[1,size(per.point.profile,2)]);
per.point(2)=psol;
varargout=dde_setupoutput('setup',funcs,per,suc,options.outputfuncs);
end
%%

