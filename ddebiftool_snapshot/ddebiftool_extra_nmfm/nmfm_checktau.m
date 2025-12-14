%% check if parameter deviations are delays
% if yes, treat as stat-dependent delay DDE
function funcs=nmfm_checktau(funcs,point,devs)
if ~funcs.tp_del
    itau=funcs.sys_tau();
    ndim=size(point.x,1);
    for i=length(devs):-1:1
        dev0(:,i)=nmfm_dev_call(devs(i),0);
    end
    pdevs=dev0(ndim+1:end,:)';
    if any(any(pdevs(:,itau)))
        %% we need delay functions
        funcs.tp_del=true;
    end
end
end
