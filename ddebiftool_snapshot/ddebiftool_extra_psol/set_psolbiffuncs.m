function extfuncs=set_psolbiffuncs(funcs,ip,rhs,method,varargin)
default={{'hjac','hdev','deviation'},0};
options=dde_set_options(default,varargin,'pass_on');
user_lhs_num=funcs.lhs_matrixfun(ip.dim);
n_ext=ip.extdim/ip.dim;
lhs_num=kron(eye(n_ext+1),user_lhs_num);
orig_xpattern=dde_funcs_xpattern(funcs,ip.dim);
xdpattern=[orig_xpattern;ones(1,size(orig_xpattern,2))*2];
xpattern=dde_join_xpattern(xdpattern,orig_xpattern);
vpattern=orig_xpattern;
for i=1:n_ext
    vpattern(1,:)=vpattern(1,:)+ip.dim;
    xpattern=dde_join_xpattern(vpattern,xpattern);
end
fun_args={'x_vectorized',true,'p_vectorized',true,...
    'lhs_matrix',lhs_num,'xpattern',xpattern};
% indices of additional delays and relations needed for extended system:
[ip.orig_ntau,deriv_order]=dde_num_delays(funcs);
if ~funcs.tp_del % constant delay
    %% set up functions of extended system
    ip.orig_tau=funcs.sys_tau();
    sys_tau=@()ip.orig_tau;
    sd_args={'sys_tau',sys_tau};
else % state-dependent delay
    %% set up functions of extended system
    sys_tau=@(it,x,p)funcs.wrap_tau(it,x(1:ip.dim,:,:),p(1,1:ip.nuserpar,:));
    sys_dirdtau=arrayfun(@(k){@(it,x,p,dx,dp)funcs.wrap_dirdtau{k}(it,...
        x(1:ip.dim,:,:),p(1,1:ip.nuserpar,:),dx(1:ip.dim,:,:),dp(1,1:ip.nuserpar,:))},...
        1:length(funcs.wrap_dirdtau));
    sd_args={'wrap_tau',sys_tau,'wrap_dirdtau',sys_dirdtau,...
        'sys_ntau',funcs.sys_ntau,'sys_tau_seq',funcs.sys_tau_seq};
end
%% required amendments of structures for extended system
sys_rhs=@(x,p)rhs(ip,funcs,0,x,p);
sys_dirderi=@(x,p,dx,dp)rhs(ip,funcs,1,x,p,dx,dp);
rhs_dirderi={'wrap_dirderi',{sys_dirderi}};
if isnumeric(options.hjac)&&options.hjac>0
    fun_dev={'hjac',options.hjac};
    rhs_dirderi={};
else
    fun_dev={};
end
extfuncs=set_funcs('wrap_rhs',sys_rhs,rhs_dirderi{:},...
    'delayed_derivs',max(deriv_order)+2,...
    sd_args{:},fun_args{:},fun_dev{:});
%% extended delays are for derivatives
embed=@(p,component,template)dde_coll_extendblanks(p,template);
extfuncs.embed=embed;
extfuncs.userfuncs=funcs;
extfuncs.usermethod=method;
extfuncs.ip=ip;
end