%% Generate right-hand side for position control problem
%
% The problem originates from H.-O. Walther, "Stable periodic motion of a
% system with state- dependent delay," Differential and Integral Equations
% 15, 923– 944 (2002).
%
% A bat or submarine wants to control its position $x(t)$ to be at position
% $x_0$ using feedback control with gain $k$. Inside the feedback loop
% there is a fixed processing delay $\tau_0$. The "current position", as
% measured is called $x_m(t)$:
%
% $$ x'(t)=k(x_0-x_m(t-\tau_0)) $$
%
% The measurement is done by sending a sound wave to the wall at position
% $0$ and measuring how long it takes for the signal to arrive back. The
% quantity $s(t)$ is the traveling time of the signal to the wall and back.
% The factor $c$ is the speed of sound such that the position estimate
% $x_m(t)$ is
%
% $$ x_m(t)=\frac{c}{2}s(t) $$
%
% The sounds true traveling time of the signal arriving at time $t$ depends
% on the sum of the position when the signal was emitted, $x(t-s(t))$ and
% the current position, $x(t)$:
%
% $$ c s(t)=x(t-s(t))+x(t) $$
%
% We generate the right-hand side symbolically, storing the output  in the
% file |sym_poscontrol.m|.
%% Load paths for symbolic r.h.s. generation
clear
base=[pwd(),'/../../'];
addpath([base,'ddebiftool']);
addpath([base,'ddebiftool_extra_symbolic']);
if dde_isoctave()
   pkg load symbolic
end
%% States
% We implement the above 3 equations as a DDAE, using the variables $x$,
% $x_m$ and $s$. We have 2 delays: $\tau_0$ and $s(t)$, such that each
% variable has 3 time instances.
%% Set number of delays and create parameter names as strings
y=sym('y',[2,3]); % symbols for states and delayed states
x=y(1,:);  % position x
s=y(2,:);  % travel time of sound s
parnames={'x0','tau0','k','c'};
cind=[parnames;num2cell(1:length(parnames))];
ip=struct(cind{:});
syms(parnames{:});     % create symbols for parameters
par=sym(parnames);
%% DDAE
% delays are in order (tau0,s) such that states are in order, e.g.,
% |x=[x(t),x(t-tau0),x(t-s(t))|. The derivatives on the left-hand side will
% be multiplied by the matrix
%
% $$ M=\left[\begin{array}{cc} 1  & 0\\ 0& 0\end{array}\right] $$
tau=[tau0;s(1)];
xm=s*c/2;
dxdt=k*(x0-xm(2));    % 
eq_s=c*s(1)-x(3)-x(1);
%% generate code
[fstr,fderivs]=dde_sym2funcs(...
    [dxdt;eq_s],... % DDAE equations
    y,...                 % variables
    par,...               % parameters
    'sd_delay',tau,...    % delays, needed only for state-dependent delays
    'filename','sym_poscontrol',...
    'sd_delay_seq',{1:2});
