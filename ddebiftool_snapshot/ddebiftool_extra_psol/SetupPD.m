%% SetupPD - Initialize continuation of period doubling bifurcation
%%
function varargout=SetupPD(funcs,branch,ind,varargin)
%% 
% Simple wrapper around SetupTorusBifurcation to have a sensible name and
% pass on bifurcation type to pdfuncs.
% See <SetupTorusBifurcation.html> for description of input and output.
%
varargout=cell(1,3);
[varargout{:}]=SetupTorusBifurcation(funcs,branch,ind,'biftype','PD','TorusInit.closest',-1,varargin{:});
end
