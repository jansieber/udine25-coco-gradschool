function [vtvres,vtvJ]=sys_cond_POEV1_norm(pt,ip,free_par_ind,varargin)
%% obtain condition that nullvector/eigenvector has length sqrt(res)
default={'res',1,'period',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
vpt=p_axpy(0,pt,[]);
vpt.profile(ip.dim+1:end,:)=pt.profile(ip.dim+1:end,:);
vpt.parameter=pt.parameter;
free_par_ind=setdiff(free_par_ind,ip.periodnullpar);
[vtv,vtvJ]=p_dot(vpt,vpt,'free_par_ind',free_par_ind,pass_on{:});
vtvJ=p_axpy(2,vtvJ,[]);
hasper=double(options.period);
vtvres=vtv+pt.parameter(ip.periodnullpar)^2+hasper*pt.period^2-options.res;
vtvJ.parameter(ip.periodnullpar)=2*pt.parameter(ip.periodnullpar);
vtvJ.period=hasper*2*pt.period;
end
