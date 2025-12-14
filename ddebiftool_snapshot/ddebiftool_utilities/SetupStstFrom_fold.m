%% Branch off at Branch point (point is either psol near BP or BP)
%%
function varargout=SetupStstFrom_fold(funcs,ststbranch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user
% * |ststbranch|: branch of |'stst'| steady state solutions or |'fold'|
% (but actually pitchfork) solutions from which one wants to branch off
% * |ind|: index in |'point'| field of |ststbranch| which is close to
% pitchfork point where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |'fold.correction'|: whether fold point should be corrected
% * |'contpar'|: index of continuation parameters if different from
% |'ststbranch'|
% * |'corpar'|: parameters left free for initial correction (if different
% from |contpars|)
%
% All other optional inputs are passed on to fields of out branch |per|
%% Outputs
%
% * (|funcs|: modified rhs functions. This output is given only
% if requested by optional argument |'outputfuncs', true|.)
% * |stst|: branch of stst solutions with desired settings and two initial
% corrected points
% * |suc|: flag indicating success
%
%%
default={'radius',0.01,'contpar',[],'corpar',[],'usercond',[],'initcond',[],...
    'fold_correction',true,'outputfuncs',false};
[foldargs,args]=dde_options_filter('SetupFold',varargin);
[options,pass_on]=dde_set_options(default,args,'pass_on');
% create branch per of periodic solutions branching off from an
% approximate BP point ind on a branch of periodic orbits (or POfold's)
if isempty(options.contpar)
    options.contpar=ststbranch.parameter.free;
end
if isempty(options.corpar)
    % assume that initial correction have to be made in all continuation
    % parameters (to accomodate step condition)
    options.corpar=options.contpar;
end
if options.fold_correction
    [foldfuncs,foldbr,suc]=SetupFold(funcs,ststbranch,ind,'outputfuncs',true,...
        'contpar',options.contpar,'corpar',options.corpar,foldargs{:},pass_on{:}); %#ok<ASGLU>
    if suc==0
        varargout=dde_setupoutput(funcs,[],suc,options.outputfuncs);
        return
    end
    fold=foldbr.point(1);
elseif isfield(ststbranch.point(ind),'v')
    fold=ststbranch.point(ind);
else
    fold=dde_fold_from_stst(ststbranch.point(ind),'funcs',funcs,pass_on{:});
end
deg_stst=dde_stst_create('point',fold);
dev=dde_stst_create('x',fold.v,'parameter',0*fold.parameter);
devnorm2=p_dot(dev,dev);
dev_stst=p_axpy(options.radius/sqrt(devnorm2),dev,deg_stst);
offststbr=replace_branch_pars(setfield(ststbranch,'point',[deg_stst,dev_stst]),...
    options.contpar,pass_on); %#ok<SFLD>
funcs=dde_add_cond('SetupStstFrom_fold',funcs,options,deg_stst);
[stst,suc]=p_correc(funcs,dev_stst,options.contpar,dev,offststbr.method.point,1, ...
    dev_stst);
if suc==0
    varargout=dde_setupoutput('setup',funcs,[],suc,options.outputfuncs);
    return;
end
offststbr.point(2)=stst;
varargout=dde_setupoutput('setup',funcs,offststbr,suc,options.outputfuncs);
end
%%

