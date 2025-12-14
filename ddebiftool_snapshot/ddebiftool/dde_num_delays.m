function [out1,out2,out3]=dde_num_delays(funcs,output)
%% return numder and type of delays in problem
if funcs.tp_del
    ntau=funcs.sys_ntau();
    ntaucalls=length(funcs.sys_tau_seq);
else
    ntau=length(funcs.sys_tau());
    ntaucalls=1;
end
max_rhs_xderiv=funcs.delayed_derivs;
if nargin==1
    output='';
end
switch output
    case {'max_rhs_xderiv','delayed_derivs','max_xderiv'}
        [out1,out2,out3]=deal(max_rhs_xderiv,ntau,ntaucalls);
    case {'ntaucalls','ncalls'}
        [out1,out2,out3]=deal(ntaucalls,ntau,max_rhs_xderiv);
    otherwise
        [out1,out2,out3]=deal(ntau,max_rhs_xderiv,ntaucalls);
end
end
