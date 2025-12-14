%% Bogdanov-Takens demo - generate right-hand side and derivatives from symbolic expressions
% 
% Demo illustrating how to branch off a transcritical Bogdanov-Takens bifurcation point
%
%% Differential equations
% 
% From: W. Jiang and Y. Yuan
% Bogdanov–Takens singularity in van der Pol’s oscillator with delayed feedback
% neural system with delay.
% Physica D: Nonlinear Phenomena, 227 (2007), pp. 149–161.
%
% $$\dot ẋ_1 = x_2,$$
%
% $$\dot x_2 = \epsilon g(x_1(t-\tau))-\epsilon(x_1^2-1)x_2-x_1,$$
%
% where
%
% g(x) = \frac{e^x-1}{c_1e^x + c_2}.

%% Add paths and load sym package if GNU Octave is used
clear
ddebiftoolpath='../../../../';
addpath(strcat(ddebiftoolpath,'ddebiftool'),...
    strcat(ddebiftoolpath,'ddebiftool_extra_symbolic'));
if dde_isoctave()
    pkg load symbolic
end
%% Create parameter names as strings and define fixed parameters
% The demo has the parameters |epsilon|, |tau| and |vartau|
parnames={'epsilon','tau','vartau'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either epsilon or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for beta, alpha and tua
par=cell2sym(parnames);  % now beta is par(1) etc
%% Define system using symbolic algebra
syms x1 x1t x2 x2t % create symbols for u1(t) u1(t-tau), u2(t), u2(t-tau)
c1 = 1/4; c2 = 1/2;
g = @(x) (exp(x) - 1)/(c1*exp(x) + c2); 
dx1_dt = tau*x2;
dx2_dt = tau*(epsilon*g(x1t) - epsilon*(x1^2 - 1)*x2 - x1);
model = [dx1_dt; dx2_dt];
vars = [x1, x1t; x2, x2t];
modelname = 'vdpo';
%% Differentiate and generate code, exporting it to sym_FHN_mf (multi-linear forms)
dde_sym2funcs(model, vars, par, 'filename',strcat('sym_',modelname,'_mf'),'directional_derivative',false); 
%% Differentiate and generate code, exporting it to sym_FHN (directional derivatives)
dde_sym2funcs(model, vars, par, 'filename',strcat('sym_',modelname),'directional_derivative',true);
