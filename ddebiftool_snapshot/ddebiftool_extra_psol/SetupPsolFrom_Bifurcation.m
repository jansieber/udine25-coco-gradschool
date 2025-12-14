%% Branch off at bifurcation point (point is either period doubling or branch point)
%%
function varargout=SetupPsolFrom_Bifurcation(funcs,psolbranch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user (or period doubling
% funcs)
% * |psolbranch|: branch of |'psol'| psol solutions from which one wants to
% branch off (could be period doubling or POfold branch)
% * |ind|: index in |'point'| field of |psolbranch| which is equal or close to
% point where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |'PeriodDoubling.correction'|: whether PeriodDoubling point should be corrected
%  alternatively |'POEV1.correction'|: whether POEV1 point should be corrected
% * |'contpar'|: index of continuation parameters if different from
% |'psolbranch'|
% * |'corpar'|: parameters left free for initial correction (if different
% from |contpars|)
%
% All other optional inputs are passed on to fields of out branch |per|
%% Outputs
%
% * (|pfuncs|: modified rhs functions. This output is given only
% if requested by optional argument |'outputfuncs', true|.)
% * |per|: branch of periodic orbits with desired settings and two initial
% corrected points
% * |suc|: flag indicating success
%
%%
default={'radius',0.01,'contpar',[],'corpar',[],'usercond',[],...
    'outputfuncs',false,'output','ini',...
    'bifname','POfold','resonance',[1,1]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[bifargs,pass_on]=dde_options_filter(['Setup',options.bifname],pass_on);
% create branch per of periodic solutions branching off from an
% approximate Period doubling point ind on a branch of periodic orbits (or Period doublings)
if isempty(options.contpar)
    options.contpar=psolbranch.parameter.free;
end
if isempty(options.corpar)
    % assume that initial correction have to be made in all continuation
    % parameters (to accomodate step condition)
    options.corpar=options.contpar;
end
if isfield(funcs,'kind')&& strcmp(funcs.kind,options.bifname)
    [biffuncs,bifpt]=deal(funcs,psolbranch.point(ind));
else
    out=cell(1,3);
    [out{:}]=feval(['Setup',options.bifname],funcs,psolbranch,ind,'outputfuncs',true,...
        'corpar',options.corpar,'output',options.output,...
        bifargs{:},pass_on{:});
    if strcmp(options.output,'sv')
        varargout=out;
        return
    end
    [biffuncs,bifbr,suc]=deal(out{:});
    if suc==0
        varargout=dde_setupoutput(funcs,[],suc,options.outputfuncs);
        return
    end
    bifpt=bifbr.point(1);
end
deg_psol=biffuncs.get_comp(bifpt,'solution');
[dev_psol,tangent]=dde_psol_from_psol(bifpt,'funcs',biffuncs,'radius',options.radius,...
    'resonance',options.resonance,'initmethod','svd',pass_on{:});
per=replace_branch_pars(setfield(psolbranch,'point',[deg_psol,dev_psol]),...
    options.contpar,pass_on); %#ok<SFLD>
funcs=dde_add_cond(['SetupPsolFrom_',options.bifname],funcs,options,dev_psol);
[psol,suc]=p_correc(funcs,dev_psol,options.contpar,tangent,per.method.point,1, ...
    dev_psol);
if suc==0
    varargout=dde_setupoutput('setup',funcs,[],suc,options.outputfuncs);
    return;
end
per.point(2)=psol;
varargout=dde_setupoutput('setup',funcs,per,suc,options.outputfuncs);
end
%%

