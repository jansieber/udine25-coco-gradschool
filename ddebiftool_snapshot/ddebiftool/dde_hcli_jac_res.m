function varargout=dde_hcli_jac_res(funcs,pt,free_par,method,varargin)
%% Jacobian and residual for connecting orbits
% wrapper distinguishing two implementations of boundary conditions
% function [J,res]=dde_hcli_jac_res(funcs,pt,free_par,method,pref,varargin)
% INPUT:
%   funcs problem functions
%	method: method parameters
%   point components
%	T period 
%	profile profile in R^(n x m*l+1)
%	t representation points in [0,1]^(m*l+1)
%	deg degree piecewise polynomial
%	par current parameter values in R^p
%	free_par free parameters numbers in N^d 
%	ph use phase condition or not (1 or 0)
%       lambda_v unstable eigenvalues of x1 (in R^s1)
%       lambda_w unstable eigenvalues of x2 (in R^s2)
%       v unstable eigenvectors of x1 (in R^n x s1)
%       w unstable eigenvectors of x2 (in R^n x s2)
%       alpha coefficients of initial functionsegment (in R^s1)
%       epsilon global coefficient of initial function segment (in R)
%       x1 steady state solution at t=-inf (in R^n)
%       x2 steady state solution at t=+inf (in R^n)
%       (optional) previous: previous solution, used in phase condition
% OUTPUT: 
%	J   jacobian in
%            R^(n*m*l+3*n+(s1+s2)*n+s1+2*s2+1 x n*m*l+3*n+(s1+s2)*n+2*s1+s2+d)
%	res residual in R^(n*m*l+3*n+(s1+s2)*n+s1+2*s2+1)
%
%   ieq structure describing the residual components
%
%
%%
varargout=cell(1,nargout);
if isfield(method,'bctype')&& strcmp(method.bctype,'subspace')
    [varargout{:}]=dde_hcli_subspace_jac_res(funcs,pt,free_par,method,varargin{:});
else
    [varargout{:}]=dde_hcli_eigval_jac_res(funcs,pt,free_par,method,varargin{:});
end
end
