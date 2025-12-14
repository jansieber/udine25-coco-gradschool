%% Branch off at Branch point (point is either period doubling or near period doubling)
%%
function varargout=SetupPsolFrom_PeriodDoubling(funcs,psolbranch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user (or period doubling
% funcs)
% * |psolbranch|: branch of |'psol'| psol solutions from which one wants to
% branch off (could be period doubling branch)
% * |ind|: index in |'point'| field of |psolbranch| which is equal or close to
% point where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |'PeriodDoubling.correction'|: whether PeriodDoubling point should be corrected
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
varargout=cell(1,nargout);
[varargout{:}]=SetupPsolFrom_Bifurcation(funcs,psolbranch,ind,...
    'fullname','PeriodDoubling','funcsname','PD','resonance',[1,2],varargin{:});
end
%%

