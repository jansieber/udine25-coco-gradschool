%% Dynamics of a inverted pendulum balancing on a cart with delayed PD feedback
%
% This demo illustrates continuation of symmetry-breaking bifurcations for
% equilibria and peridoic orbits.
%
% This file generates the right-hand side generated using symbolic toolbox
% in |gen_sym_balancing|.
%
% Differential equations are 
% 
% $$x''(t)=\sin x(t)-\cos x(t)[ax(t-\tau)+bx'(t-\tau)]$$
%
% after non-dimensionalization.
%
% Reference:
% 
% J. Sieber, B. Krauskopf, Bifurcation analysis of an inverted pendulum
% with delayed feedback control near a triple-zero eigenvalue. Nonlinearity
% 17 (1), pp. 85-104, 2004.
%
%%
clear
addpath([pwd(),'/../../ddebiftool']);
addpath([pwd(),'/../../ddebiftool_extra_symbolic']);
if dde_isoctave()
    pkg load symbolic
end
%% Set number of delays and create parameter names as strings
syms x x_tau v v_tau % symbols for states and delayed states
parnames={'a','b','tau'};
cind=[parnames;num2cell(1:length(parnames))];
ind=struct(cind{:});
syms(parnames{:});     % create symbols for parameters
par=sym(parnames);
%% DDE
dxdt=v;
dvdt=sin(x)-cos(x)*(a*x_tau+b*v_tau);
%% generate code
[fstr,fderivs]=dde_sym2funcs(...
    [dxdt;dvdt],...
    [x,x_tau;v,v_tau],...
    par,...
    'directional_derivative',true,'filename','sym_balancing');
