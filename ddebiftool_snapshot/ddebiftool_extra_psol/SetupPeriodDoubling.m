%% SetupPeriodDoubling - Initialize continuation of period doubling bifurcation
%%
function varargout=SetupPeriodDoubling(funcs,branch,ind,varargin)
%% 
% Simple wrapper around SetupPD.
% See <SetupTorusBifurcation.html> for description of input and output.
%
varargout=cell(1,nargout);
[varargout{:}]=SetupPD(funcs,branch,ind,varargin{:});
end
