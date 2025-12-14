%% Generaluzed-Hopf demo - generate right-hand side and derivatives from symbolic expressions
%
% <html>
% $Id$
% </html>
% 
% From:
% Zhen, Bin (PRC-TONG-AEM); Xu, Jian [Xu, Jian2] (PRC-TONG-AEM)
% Bautin bifurcation analysis for synchronous solution of a coupled FHN
% neural system with delay.
% Commun. Nonlinear Sci. Numer. Simul. 15 (2010), no. 2, 442–458.
%
%% Differential equations
%
% $$u_1'=-\frac{u_1^3+(c+\alpha)u_1^2+d u_1-u_2+2\beta*tanh(u_1(t-\tau))$$
%
% $$u_2'=\epsilon(u_1-b u_2)$$
%% Add paths and load sym package if GNU Octave is used
clear
ddebiftoolpath='../../../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_symbolic'));
if dde_isoctave()
    pkg load symbolic
end
%% Create parameter names as strings and define fixed parameters
% The demo has the parameters |beta|, |alpha| and |tau|
parnames={'Q','E','tau'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either beta or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for beta, alpha and tua
par=cell2sym(parnames);  % now beta is par(1) etc
%% Define system using symbolic algebra
syms u1 u1t u2 u2t % create symbols for u1(t) u1(t-tau), u2(t), u2(t-tau)
alpha = @(u1) 1/(1+exp(-4*u1)) - 1/2; 
du1_dt = -u1 + 2.6*alpha(u1t) - Q*u2t + E;
du2_dt = -u2 + alpha(u1t);
model = [du1_dt; du2_dt];
vars = [u1,u1t;u2,u2t];
modelname = 'neural_network';
%% Differentiate and generate code, exporting it to sym_FHN_mf (multi-linear forms)
dde_sym2funcs(model, vars, par, 'filename',strcat('sym_',modelname,'_mf'),'directional_derivative',false); 
%% Differentiate and generate code, exporting it to sym_FHN (directional derivatives)
dde_sym2funcs(model, vars, par, 'filename',strcat('sym_',modelname),'directional_derivative',true);
