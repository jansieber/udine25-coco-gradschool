%% construct eigenvalue problem for Floquet multipliers using sparse matrices
% for case with only one-sided delays. This is based on eigs and permits
% optional argument 'closest' for finding Floquet multipliers close to a
% particular value.
%
% $Id: dde_psol_mult_app_sparse.m 352 2019-06-19 16:03:03Z jansieber $
%%
function [s,ef]=dde_psol_mult_app_sparse(Margs,varargin)
default={'geteigenfuncs',false,...
    'closest',[],'max_number_of_eigenvalues',20,'method',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on','method');
if isnumeric(Margs{1})
    dim=size(Margs{1},1);
    M_is_fun=false;
else
    dim=Margs{2};
    M_is_fun=true;
end
n_ev=min(options.max_number_of_eigenvalues,dim);
closest_eigs=num2cell(options.closest);
if isempty(closest_eigs) || M_is_fun
    closest_eigs={'lm'};
end
if ~options.geteigenfuncs
    s=eigs(Margs{:},n_ev,closest_eigs{:});
    ef=[];
else
    [ef,s]=eigs(Margs{:},n_ev,closest_eigs{:});
    s=diag(s);
end
if M_is_fun && ~isempty(options.closest)
    s=options.closest+1./s;
end
end
