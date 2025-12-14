function binv = nmfm_binv(funcs, point, lambda, q, p, zeta, kappa)
%% Compute a bordered inverse
% INPUT:
%   point: stst type point (local bifurcation)
%	lambda: eigenvalue of the characteristic matrix Delta
%   q: normalized null vector of Delta
%   p: normalized null vector of Delta^T
%   zeta: n-vector
%   kappa: scalar
% OUTPUT:
%	binv: solution history function to the bordered inverse problem
%
% $Id: nmfm_binv.m 314 2019-01-24 14:28:23Z mmbosschaert $
%
%% characteristic matrix
chmatfun=@(deg)ch_matrix(funcs,point.x,point.parameter,lambda,'deri',deg);
[Delta0,Delta1,Delta2]=deal(chmatfun(0),chmatfun(1),chmatfun(2));
%% solution of bordered problem
% buggy? should A=[Delta0,Delta1*q; p*Delta1,0], because A as below can be
% singular if p*q=0.
A = [Delta0, q; p, 0];
B = [zeta + kappa*Delta1*q; 0];
X = A\B;
xi = X(1:end-1);
gamma = -p*Delta1*xi + (1/2)*kappa*p*Delta2*q;
p0=zeros(size(point.parameter(:)));
binv=nmfm_dev_fun( [xi+gamma*q,-kappa*q;p0,p0], 'lambda',lambda*[1,1],'t',[0,1]);
end

