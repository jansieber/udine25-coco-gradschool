function branch=df_brnch(funcs,free_par,kind,varargin)
%% set up empty branch with default values
% function branch=df_brnch(funcs,free_par,kind,flag_newhheur)
% INPUT:
%       funcs problem functions
%       free_par free parameter list
%       kind type of solution point
%       flag_newhheur (optional, default: 2) use new Chebyshev approximation
%                     Only used if kind==stst, cf. df_mthod
% OUTPUT:
%	branch empty branch with default values

% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
% Update on 05/03/2007 ("flag_newhheur" <=> (imag(method.stability.lms_parameter_rho)~=0) )
%
% $Id: df_brnch.m 296 2018-09-24 21:36:56Z jansieber $
%
%%


branch.method=df_mthod(kind,varargin{:});

branch.parameter.free=free_par;
branch.parameter.min_bound=[];
branch.parameter.max_bound=[];
branch.parameter.max_step=[];

if isempty(funcs) || ~isstruct(funcs)
    return
end
branch=br_bound_delay(funcs,branch);
branch.point=[];
end
