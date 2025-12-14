%% call derivatives of delays according to sys_tau_seq
function dtaus=dde_dtau_dirderi(funcs,x,p,dx,dp,varargin)
default={'order',1,'repeat',true};
options=dde_set_options(default,varargin,'pass_on');
x=x(:,:,:,1);
dx=dx(:,:,:,1);
dtaus=dde_all_dirdtau(funcs,true,options.order,x,p,dx,dp);
if options.repeat
    [ntau,deriv_order]=dde_num_delays(funcs);
    dtaus=reshape(repmat(dtaus,[deriv_order+1,1]),...
        (ntau+1)*(deriv_order+1),[]);
end
end
