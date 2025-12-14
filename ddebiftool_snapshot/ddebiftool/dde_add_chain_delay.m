function ind=dde_add_chain_delay(ind,name,g,varargin)
%% add nested distributed delay by adding equation of the form
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
% idist are index of d in variable array, ibd is
% index of tau in variable array, ginp is kernel function g(s,x,p), ntau is
% the delay number after which the new delays (nodes in the integral)
% should be appended (ntau+(1:length(grid.t))
%
% optional argument 'par' passes index of parameters in argument p for g
% (default same as for DDE rhs), 'gt' or 'g_of_t_x flags that first argument t
% of g is present (default true)
%% wrap g, find its image dimension
[g,ng]=dde_sym_dist_delay(g,varargin{:});
default={{'igx','x_args','x_arguments'},NaN,...
    {'igp','p_args','p_arguments','ipar'},zeros(1,0),{'ig0','igini'},NaN,...
    {'ig0type','iginitype'},'0'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% add name to list of names
num=length(fieldnames(ind.histories.names));
ind.histories.names.(name)=num;
ind.histories.fcn{num}=@(varargin)dir_deriv(g{1},g(2:end),[1,1,2],...
    varargin{:},pass_on{:},'nf',ng);
iname=ind.histories.names.(name);
ind.histories.xvals{iname+1}{ind.histories.int}=ind.histories.lastxarg+(1:ng);
ind.histories.xvals{iname+1}{ind.histories.id}= ind.histories.lastxarg+ng+(1:ng);
ind.histories.xinit{iname}=ind.histories.lastinit+(1:ng);
ind.histories.lastxarg=ind.histories.lastxarg+2*ng;
ind.histories.lastinit=ind.histories.lastinit+ng;
%% set up starting point for integral
ind.args.init.global{num}=options.ig0(:).';
if iscell(options.ig0type)
    ind.args.init.type{num}=options.ig0type(:)';
else
    ind.args.init.type{num}=repmat({options.ig0type},1,length(options.ig0));
end
%% collect x arguments of g
% format igx argument
ind.histories.xinputs{num}=int_delay_chain_inputs(ind.histories,options.igx);
%% collect p arguments for g
ind.args.igp.global{num}=options.igp(:).';
ind.args.igp.type{num}=repmat({'parameter'},1,length(options.igp));
ind.msh.wJ{num}=kron(speye(ng),ind.msh.Jibase);
end
