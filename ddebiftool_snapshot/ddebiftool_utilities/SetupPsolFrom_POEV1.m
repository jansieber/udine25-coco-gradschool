%% Branch off at Branch point (point is either psol near BP or BP)
%%
function varargout=SetupPsolFrom_POEV1(funcs,psolbranch,ind,varargin)
%% Inputs
% 
% * |funcs|: structure with functions provided by user
% * |psolbranch|: branch of |'psol'| psol solutions from which one wants to
% branch off
% * |ind|: index in |'point'| field of |psolbranch| which is close to
% point where we want to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |'POfold.correction'|: whether POfold point should be corrected
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
varargout=cell(1,3);
[varargout{:}]=SetupPsolFrom_Bifurcation(funcs,psolbranch,ind,...
    'fullname','POEV1','funcsname','POfold','resonance',[1,1],varargin{:});
end
%%

