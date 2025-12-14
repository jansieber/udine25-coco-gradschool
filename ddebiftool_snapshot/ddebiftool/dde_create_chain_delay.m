function ind=dde_create_chain_delay(bound,ix,varargin)
%% prepare nested distributed delay for adding equations of the form
% h_{j+1}(t-s) = \int_0^s g_j(s,x_j(t-s),par) ds for j=0..k
% 0 = \int_0^(tau*t(i)) h_{j_i}(s) ds - d_i 
%
% where h0=x0=x and xj contains arbitrary combination of h_i with i<j
%
% 0 = int_0^inf g(s,x(t-s),par) exp(-s/tau) (s/tau)^alpha ds - d 
% or equivalently
% 0 = int_0^inf tau g(tau s,x(t-tau s),par) exp(-s) (s)^alpha ds - d 
% (generalized) Laguerre approximation
%
% ibd is index of tau in variable array or parameter array, dependingon
% option bound_is_par. Options |nint| and |degree| fix discretization.
% Option |'infinite_delay'| treats delay as Laguerre case:
% \int_0^\infty g(s,x(t-s),p)\exp(-\tau s)ds
default={{'int','nint','tcoarse','grid'},4,'degree',3,...
    {'bdtype','bound','boundtype'},'parameter',...
    {'laguerre','infinite_delay'},false,'laguerre_alpha',0};
options=dde_set_options(default,varargin,'pass_on');
%% define mesh (non-adaptive)
if options.laguerre
    type='laguerre';
    nint=options.laguerre_alpha;
else
    type='cheb';
    nint=options.int;
end
ind.msh=dde_dist_integrator(options.degree,nint,type);
%% initialize argument type definitions
ind.args=struct(...
    'ibd', arg_create('type',{{options.bdtype}},...
                      'global',{bound}),...
    'igx', arg_create('type',{repmat({'x'},1,length(ix))},...
                      'global',{ix(:).'}),...
    'igp', arg_create(),...
    'ival',arg_create(),...
    'init',arg_create('replace',{'0',@(nx,nvec)zeros(nx,nvec)}),...
    't',   arg_create('replace',{'NaN',@(nx,nvec)NaN(nx,nvec)}));
%% initialize history segments
ind.histories.fcn={};
ind.histories.int=1;
ind.histories.id= 2;
ind.histories.names=struct('x',0);
ind.histories.xvals={{[],1:length(ix)}};
ind.histories.lastxarg=length(ix);
ind.histories.lastinit=0;
ind.histories.xinputs={};
ind.histories.xinit={};
end
%%
function arg=arg_create(varargin)
default={'global',{},'nval',0,'type',{},'local',struct(),'replace',{}};
arg=dde_set_options(default,varargin);
end