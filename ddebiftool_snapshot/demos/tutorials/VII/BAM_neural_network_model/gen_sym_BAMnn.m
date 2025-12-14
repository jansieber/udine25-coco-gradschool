%% Transcritical Bogdanov-Takens demo - generate right-hand side and derivatives from symbolic expressions
%
% <html>
% $Id$
% </html>
% 
% From:
% From: Dong, Tao and Liao, Xiaofeng
% Bogdanov-Takens bifurcation in a tri-neuron BAM neural network model with multiple delays
% Nonlinear Dynamics.(2013), no. 3, 583-595.
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
% The demo has the parameters |alpha1|, |alpha2| and |tau|
parnames={'alpha1','alpha2','tau0','tau'};
%% Create symbols for parameters, states and delays states
% The array |par| is the array of symbols in the same order as parnames.
% Due to the following two lines we may, for example, use either beta or
% par(1) to refer to the delay.
syms(parnames{:});       % create symbols for beta, alpha and tua
par=cell2sym(parnames);  % now beta is par(1) etc
%% Define system using symbolic algebra
syms u1 u1t u2 u2t u3 u3t% create symbols for u1(t) u1(t-tau), u2(t), u2(t-tau)
mu1 = 0.1; mu2 = 0.3; mu3 = 0.2; c12 = 1; c13 = 1;
f1 = @(x) tanh(x) + 0.1*x^2;
f2 = @(x) tanh(x);
f3 = @(x) tanh(x);
c210 =  (mu2^2*(mu1 + mu3 + mu1*mu3*tau0))/(c12*(mu2 - mu3))
c310 = -(mu3^2*(mu1 + mu2 + mu1*mu2*tau0))/(c13*(mu2 - mu3))
% c210 = 0.36;
% c310 = -0.22;
du1_dt = -mu1*u1 + (c210 + alpha1)*f1(u2t) + (c310 + alpha2)*f1(u3t);
du2_dt = -mu2*u2 + c12*f2(u1);
du3_dt = -mu3*u3 + c13*f3(u1);
modelname = 'BAMnn';
model = [du1_dt;du2_dt;du3_dt];
vars  = [u1,u1t;u2,u2t;u3,u3t];
%% Differentiate and generate code, exporting it to sym_FHN_mf (multi-linear forms)
dde_sym2funcs(model, vars, par, 'filename',strcat('sym_',modelname,'_mf'),'directional_derivative',false); 
%% Differentiate and generate code, exporting it to sym_FHN (directional derivatives)
dde_sym2funcs(model, vars, par, 'filename',strcat('sym_',modelname),'directional_derivative',true);
