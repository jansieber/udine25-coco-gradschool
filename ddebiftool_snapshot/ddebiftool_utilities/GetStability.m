%% Number of unstable eigenvalues/Floquet multipliers along branch
%%
function [nunst,dom,triv_defect,points]=GetStability(branch,varargin)
%% Inputs
%
% * |branch|: branch along which number of unstable eigenvalues/Floquet
% multipliers is determined
%
%  Important name-value pair inputs
% 
% * |'exclude_trivial'| (logical): exclude trivial eigenvalues (what they are
%                            depends on the type of points)
% * |'locate_trivial'| (function): e.g. @(p)[1;-1] for removing 1 and -1
%                            for period doubling orbits (overwrites standard location)
% * |'funcs'|:   set of functions for computing branch, needed if stability 
% information is not yet present and needs to be calculated
% * |'method'|: which method is used for computation (if not present, a
% standard method is generated with |df_mthod|)
% * |'recompute'| (default |false|):  force recomputation of stability even
% if it is present.
%
%% Outputs
% * |nunst| (vector of integers): number of unstable eigenvalues
% * |dom| (vector of real/complex): dominant eigenvalue (closest unstable to
%                       bifurcation if exists, otherwise, closest stable)
% * |triv_defect| (vector of real): distance of eigenvalue approximating
%                               trivial eigenvalue from its true value
% * |points| (struct array): array of points (same as in input branch but
% with stability info if not present on input)
%
% $Id: GetStability.m 369 2019-08-27 00:07:02Z jansieber $
%
%% process options
defaults={'funcs',[],'exclude_trivial',false};
[options,pass_on]=dde_set_options(defaults,varargin,'pass_on');
%% check if stability has been computed already & compute if necessary
if ~isfield(branch,'point')
    points=branch;
    branch=df_brnch(options.funcs,[],points(1).kind);
    branch.point=points;
    branch.method.stability=dde_set_options(branch.method.stability,pass_on,'pass_on');
end
[branch,nunst,dom,triv_defect]=br_stabl(options.funcs,branch,...
    'exclude_trivial',options.exclude_trivial,pass_on{:});
points=branch.point;
end
