%% Definition of user functions
%
% <html>
% $Id: demo1_funcs.m 367 2019-07-14 21:56:46Z jansieber $
% </html>
%
% (Please load |ddebiftool| into path first, see <demo1.html>.)
% To define a system, the user should provide Matlab functions defining the
% right-hand side (in the example called |neuron_sys_rhs|) and the indices
% for the delays in the parameter vector (called |neuron_tau| in our
% example).
%% Right-hand side
% A function defining the right-hand side $f$:
% 
%   function y=sys_rhs(x,p)
% 
% This function has two arguments, |x| $\in R^{n_x\times (n_\tau+1)}$, which contains 
% the state variable(s) at the present and in the past,
% |x| $=[x(t), x(t-\tau_1), \ldots, x(t-\tau_{n_\tau})]$,
% |p| $\in R^{1\times n_p}$ which contains the parameters, |p|.
%
% For the example, this is ($f$ is called |neuron_sys_rhs| in our example)
%
% This simple demo manually numbers the parameters in the array as
% |par(1)...par(7)|. See <demo1_simple.html> for a safer way to keep
% track of parameter ordering in the parameter array.
neuron_sys_rhs_nonvec=@(x,p)[...
    -p(1)*x(1,1)+p(2)*tanh(x(1,4))+p(3)*tanh(x(2,3));....
    -p(1)*x(2,1)+p(2)*tanh(x(2,4))+p(4)*tanh(x(1,2))];
%% Vectorization
% Computations will be faster (especially for periodic orbits) if the
% function can be called in many x values and parameters at once.
neuron_sys_rhs=@(x,p)[...
    -p(1,1,:).*x(1,1,:)+p(1,2,:).*tanh(x(1,4,:))+p(1,3,:).*tanh(x(2,3,:));....
    -p(1,1,:).*x(2,1,:)+p(1,2,:).*tanh(x(2,4,:))+p(1,4,:).*tanh(x(1,2,:))];

%% Delays
% The delays $\tau_i$, $i=1\ldots,m$ are considered to be part of the
% parameters ($\tau_i=\eta_{j(i)}$, $i=1,\ldots,m$). This is natural since
% the stability of steady solutions and the position and stability of
% periodic solutions depend  on the values of the delays. Furthermore
% delays can occur both as a 'physical' parameter and as delay, as in
% $\dot{x}=\tau x(t-\tau)$. From these inputs the right-hand side $f$ is
% evaluated at time $t$. For equations with constant delays DDE-Biftool
% determines which parameters are delays by calling an argument-less
% function of the type
%
%   function d=sys_tau()
%
% In the example we order the parameters as |par| $=[\kappa,
% \beta, a_{12}, a_{21},\tau_1,\tau_2, \tau_s]$. Thus, (giving it the name
% |neuron_tau|):
neuron_tau=@()[5,6,7];
ind_a21=4;  % used later for continuation
ind_taus=7; % used later for continuation

%% Deivatives of user-provided functions
% Optionally (recommended) the user may also specify the directional
% derivatives of the user-defined functions with respect to states, delayed
% states and parameters. For constant delays only the derivatives of $f$
% are required. They should be provided as a cell array of functions of the form
%
%   function df=sys_dirderi(xx,par,dxx,dpar)
%
% providing the derivatives of order k in directions dxx, dpar (dxx and
% dpar have the same shape as xx and par).
%
%% Definition of structure |funcs|
% Similar to standard functions such as |ode45| DDE-Biftool's routines have
% an argument that defines the right-hand side. Since DDE-Biftool needs
% several user-defined functions (sys_rhs, sys_tau, optionally sys_deri,
% sys_cond etc) these functions are collected in a structure |funcs|. This
% structure |funcs| is best set up by calling the DDE-Biftool routine
% |set_funcs| with a sequence of name-value pairs. Each name-value pair
% corresponds to a field in the structure. Fields that are not listed as
% arguments of set_funcs get replaced by a default if possible.
%
% Possible argument names are:
% 
% * |'sys_rhs'| (default |sys_rhs| if file |sys_rhs.m| present in folder):
%    right-hand side |sys_tau|
% * |'sys_tau'| (default |@()[]|): function defining delays
% * |'sys_dirderi'| (default |[]|): function defining directional
%     derivatives of |sys_rhs|. This is an alternative to providing
%     |sys_deri|. |sys_dirderi| may be a cell array, where entry
%     sys_dirderi{1}(xx,par,dx,dpar) is the first derivative in direction
%     (dx,dpar). Orders 1 and 2 are used for standdard bifurcation
%     analysis, higher orders (up to 5, only with dp=0) for normal form
%     computations
% * |'sys_ntau'| (default 0, only needed for state-dependent delays) number
%    of delays
% * |'sys_cond'| (default |@dummy_cond|) function or structure providing extra conditions
% * |'sys_dirdtau'| (default |[]|, only needed for state-dependent delays):
%      function defining directional derivatives of |sys_tau|. 
% * |'x_vectorized'| (logical, default false) set to true if |sys_rhs|,
%     |sys_dirderi|, (|sys_tau| and |sys_dirdtau| for SD-DDEs) accept an argument
%    |xx| with three dimensions. For periodic-orbit computations and
%    Jacobians the function will be called with many arguments
%     simultaneously if |x_vectorized| is true.
% * |'p_vectorized'| (logical, default false) set to true if |sys_rhs|,
%     |sys_dirderi|, (|sys_tau| and |sys_dirdtau| for SD-DDEs) accept an argument
%    |par| with three dimensions (1 x np x nvec). For periodic-orbit
%    computations and Jacobians the functions will be called with many
%    argumentss simultaneously if |x_vectorized| is true.
%
% Other fields will be set internally: |tp_del| (true if delays are state-dependent),
% |sys_dirderi_provided| (true if user has provided |sys_dirderi|) and
% |sys_dirdtau_provided| (true if user has provided |sys_dirdtau|).
% |wrap_rhs|, |wrap_...| are wrapped versions of the user-provided
% functions that are guaranteed to be callable in vectorized fashion.
% |funcs.drhs_dir(ord,x,p,dx,dp)| returns arbitrary-order directional
% derivatives, |funcs.drhs_mf(x,p,dx1,dp1,...)| returns arbitrary-order
% mixed derivatives and permits expansions such a
% |funcs.drhs_mf(x,p,{1,'I'},0)| to obtain Jacobians w.r.t x.
%% Numerical finite differences for derivatives
% Instead one my rely on numerical finite differences for the derivatives.
fnumeric=set_funcs(...
    'sys_rhs',neuron_sys_rhs,...
    'sys_tau',neuron_tau,...
    'x_vectorized',true,'p_vectorized',true) %#ok<NOPTS>
%% Directional derivatives
% One may define directional derivatives for more robust bifurcation and
% normal form computations:
dtanh=@(x)(1-tanh(x).^2);
ddtanh=@(x)2*(tanh(x).^2-1).*tanh(x);
neuron_sys_dirderi{1}=@(x,p,dx,dp)[...
    -dp(1,1,:).*x(1,1,:)-p(1,1,:).*dx(1,1,:)+...
     dp(1,2,:).*tanh(x(1,4,:))+p(1,2,:).*dtanh(x(1,4,:)).*dx(1,4,:)+...
     dp(1,3,:).*tanh(x(2,3,:))+p(1,3,:).*dtanh(x(2,3,:)).*dx(2,3,:);...
    -dp(1,1,:).*x(2,1,:)-p(1,1,:).*dx(2,1,:)+...
     dp(1,2,:).*tanh(x(2,4,:))+p(1,2,:).*dtanh(x(2,4,:)).*dx(2,4,:)+...
     dp(1,4,:).*tanh(x(1,2,:))+p(1,4,:).*dtanh(x(1,2,:)).*dx(1,2,:)];
neuron_sys_dirderi{2}=@(x,p,dx,dp)[...
    -2*dp(1,1,:).*dx(1,1,:)+...
     2*dp(1,2,:).*dtanh(x(1,4,:)).*dx(1,4,:)+p(1,2,:).*ddtanh(x(1,4,:)).*dx(1,4,:).^2+...
     2*dp(1,3,:).*dtanh(x(2,3,:)).*dx(2,3,:)+p(1,3,:).*ddtanh(x(2,3,:)).*dx(2,3,:).^2; ...
    -2*dp(1,1,:).*dx(2,1,:)+...
     2*dp(1,2,:).*dtanh(x(2,4,:)).*dx(2,4,:)+p(1,2,:).*ddtanh(x(2,4,:)).*dx(2,4,:).^2+...
     2*dp(1,4,:).*dtanh(x(1,2,:)).*dx(1,2,:)+p(1,4,:)*ddtanh(x(1,2,:))*dx(1,2,:).^2];
fdirectional=set_funcs(...
    'sys_rhs',neuron_sys_rhs,...
    'sys_tau',neuron_tau,...
    'sys_dirderi',neuron_sys_dirderi,...
    'x_vectorized',true,'p_vectorized',true) %#ok<NOPTS>
%% Derivatives generated by symbolic toolbox
% If the symbolic toolbox (7.0 and greater) is available, one may generate
% derivatives and right-hand side automatically (see script
% <gen_sym_demo1.html>). Gnu Octave's package symbolic may also work. The
% output of the symbolic codee generation needs to be wrapped. The function
% |set_symfuncs| provides this (it calls |set_funcs| internally).
fsymbolic=set_symfuncs(@sym_neuron,'sys_tau',neuron_tau);
%% Select one of the ways to define the problem for the demo
funcs=fsymbolic;
%% Save and continue to continuation and stability of steady states <demo1_stst.html>
save('demo1_funcs_results.mat');
