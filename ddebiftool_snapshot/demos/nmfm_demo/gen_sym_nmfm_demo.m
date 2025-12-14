%% Generate right-hand sides and derivatives for nmfm_demo using symbolic toolbox
%% The system
% This system represents a non-dimensionalized model of two interacting
% layers of neurons.
% 
% $$ \dot{x}_1(t) = - x_1(t) - a g(bx_1(t-\tau_1)) + cg (dx_2(t-\tau_2)), $$
%
% $$ \dot{x}_2(t) = - x_2(t) - a g(bx_2(t-\tau_1)) + cg (dx_1(t-\tau_2)), $$
%
% where $g:\bf{R} \rightarrow \bf{R}$ is of the sigmoidal form
%
% $$ g(z) = \left[\tanh(z-1) + \tanh(1)\right]\cosh(1)^2. $$
% 
% The variables $x_1(t)$ and $x_2(t)$ represent the population-averaged
% neural activity at time $t$ in layers one and two, respectively.
% The parameter $a > 0$ is a measure of the strength of inhibitory feedback,
% while $c > 0$ measures the strength of the excitatory effect of one layer
% on the other. The parameters $b > 0$ and $d > 0$ are saturation rates and
% the delays $\tau_{1,2}$ represent time lags in the inhibitory feedback
% loop and excitatory inter-layer connection. Note that the system is
% symmetric with respect to interchanging the labels $1$ and $2$, so
%%
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool'],...
    [base,'ddebiftool_extra_symbolic']);
if dde_isoctave()
    pkg load symbolic
end
parnames={'a','b','c','d','tau1','tau2'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either tau or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for parameter names
par=cell2sym(parnames);  % now a is par(1) etc
x=sym('x',[2,1]);        % current state
xt=sym('xt',[2,2]);      % delayed state
sigm=@(z)(tanh(z-1)+tanh(sym(1)))*cosh(sym(1))^2;
f=@(i,x,xt)-x(i)-a*sigm(b*xt(i,1))+c*sigm(d*xt(3-i,2));
dx=[f(1,x,xt);f(2,x,xt)];
%% Differentiate and generate code, exporting it to sym_minimal_demo
[fstr,derivs]=dde_sym2funcs(...
    dx,...             % n x 1 array of derivative symbolic expressions
    [x,xt],...         % n x (ntau+1) array of symbols for states (current & delayed)
    par,...            % 1 x np (or np x 1) array of symbols used for parameters
    'filename','sym_nmfm_demo',... % optional argument specifying output file
    'directional_derivative',true);
