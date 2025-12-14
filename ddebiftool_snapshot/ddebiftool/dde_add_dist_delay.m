function tab=dde_add_dist_delay(tab,g,varargin)
%% add distributed delay by adding equation of the form
% 0 = \int_0^tau g(s,x(t-s),par) ds - d 
% or equivalently
% 0 = \int_0^1 tau g(tau s,x(t-tau s),par) ds - d 
%
% or
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
tab=dde_funtab(tab);
[g,ng,gname]=dde_sym_dist_delay(g,varargin{:}); %#ok<ASGLU> % preprocess g from symbolic
default={'name',gname,...
    {'int','nint','tcoarse','grid'},4,'degree',3,'t_quad',[],'w_quad',[],...
    'delayfun',@int_delay_rhs,...
    {'igx','ix','x','kernelargs','ikargs'},NaN,...
    {'igp','ipar','par','parameter','kernelpar'},unique(cat(2,tab.fun_ref.par)),...
    {'ibd','bound'},NaN,...
    {'ival','value'},NaN,...
    {'bound_is_par','bound_is_parameter'},true,...
    'infinite_delay',false,'laguerre_alpha',0};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% determine integration weights and nodes
grid = loc_int_wts(options);
%% igx,igp variables are for arguments x,p of g(s,x,p), 
% ibd is negative of lower integral boundary, 
% tmsh are numbers of delays for integration nodes
gwrap=@(varargin)dir_deriv(g{1},g(2:end),[1,1,2],varargin{:},pass_on{:},...
    'nf',length(options.ival));
ind=struct('tau',1:length(grid.t),'igx',1:length(options.igx),...
    'igp',1:length(options.igp),...
    'val',length(options.igx)+(1:length(options.ival)));
if options.bound_is_par
    [ibdpar,ibdx]=deal(options.ibd,[]);
    ind.bd=length(options.igp)+1;
else
    [ibdpar,ibdx]=deal([],options.ibd);
    ind.bd=length(options.igx)+length(options.ival)+1;
end
ind.ntau=length(grid.t);
bd=@(x,p)bdparchoice(x,p,options.bound_is_par,ind);
rhs=@(x,p)options.delayfun(gwrap,grid,0,...
    x(ind.igx,1+ind.tau,:),p(1,ind.igp,:),bd(x,p),x(ind.val,1,:));
drhs=@(ord,x,p,dx,dp)options.delayfun(gwrap,grid,ord,...
     x(ind.igx,1+ind.tau,:), p(1,ind.igp,:),bd( x, p), x(ind.val,1,:),...
    dx(ind.igx,1+ind.tau,:),dp(1,ind.igp,:),bd(dx,dp),dx(ind.val,1,:));
drhsc=arrayfun(@(i){@(x,p,dx,dp)drhs(i,x,p,dx,dp)},1:length(g)-1);
tau =    @(it,x,p)      dist_delay_tau(ind.tau(it),bd(x,p),grid,0);
dtau=@(ord,it,x,p,dx,dp)dist_delay_tau(ind.tau(it),bd(dx,dp),grid,ord);
dtauc=arrayfun(@(ord){@(it,x,p,dx,dp)dtau(ord,it,x,p,dx,dp)},1:length(g)-1);
%% set xpattern
[ix,it]=ndgrid(ind.igx(:)',ind.tau+1);
xpattern=[ix(:),it(:);...
    ind.val(:),ones(length(ind.val),1);...
    reshape(ibdx,[],1), ones(length(ibdx),1)]';
%% generate funcs structure
xdep=[options.igx(:);options.ival(:);ibdx(:)];
nxdep=length(xdep);
funcs=set_funcs('wrap_rhs',rhs,...
    'wrap_dirderi',drhsc,...
    'wrap_tau',tau,...
    'wrap_dirdtau',dtauc,...
    'sys_ntau',@()ind.ntau,...
    'sys_tau_seq',{ind.tau},...
    'x_vectorized',true,'p_vectorized',true,...
    'xpattern',xpattern,...
    'lhs_matrix',zeros(length(ind.val),nxdep));
%% append new funcs to tab
tab=dde_add_funcs(tab,options.name,funcs,'par_old',[options.igp(:);ibdpar(:)],...
    'x_old',xdep,'rhstau_args',ind.ntau,...
    'nf',length(options.ival),pass_on{:});
end
%%
function bdval=bdparchoice(x,p,bound_is_par,ind)
if ~bound_is_par
    bdval=x(ind.bd,1,:);
else
    bdval=p(1,ind.bd,:);
end
end
%%
function grid = loc_int_wts(options)
if options.infinite_delay
    [grid.t,grid.w]=lagpts(options.degree,options.laguerre_alpha);
    return
end
if isempty(options.t_quad)
    [options.t_quad,options.w_quad]=legpts(options.degree,[0,1]);
end
if length(options.int)==1 && options.int==floor(options.int)
    tcoarse=linspace(0,1,options.int+1);
else
    tcoarse=options.int;
end
grid.w=dde_coll_meshfill(tcoarse,1,'grid',options.w_quad,'acc',false);
grid.t=dde_coll_meshfill(tcoarse,1,'grid',options.t_quad);
grid.degree=options.degree-1;
grid.tcoarse=tcoarse;
end
%%
function tau=dist_delay_tau(itau,bd,msh,order)
nvec=size(bd,3);
bd=reshape(bd,1,nvec);
if order<=1
    tau=reshape(msh.t(itau),[],1)*bd;
else
    tau=zeros(length(itau),nvec);
end
end
