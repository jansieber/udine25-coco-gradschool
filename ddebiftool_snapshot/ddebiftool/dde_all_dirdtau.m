function J=dde_all_dirdtau(funcs,inc0,order,x,p,dx,dp)
if ~funcs.tp_del
    dpsize=size(dp);
    vecdim=[dpsize(3:end),1];
    if order==1
        dp=reshape(dp,size(dp,2),prod(vecdim));
        J=dp(funcs.sys_tau(),:);
    else
        nvec=prod(vecdim);
        J=zeros(length(funcs.sys_tau()),nvec);
    end
    if inc0
        J=[zeros(1,size(J,2));J];
    end
else
    xsize=size(x);
    dxsize=size(dx);
    vecdim=[dxsize(3:end),1];
    psize=size(p);
    x=reshape(x,xsize(1),xsize(2),[]);
    p=reshape(p,psize(1),psize(2),[]);
    taucalls=[{1},cellfun(@(i){i+1},funcs.sys_tau_seq)];
    ncalls=length(taucalls);
    for nc_i=ncalls:-1:2
        t_i=taucalls{nc_i};
        mt_i=min(t_i)-1;
        J(t_i,:)=funcs.wrap_dirdtau{order}(t_i-1,...
            x(:,1:mt_i,:),p,dx(:,1:mt_i,:),dp);
    end
    if ~inc0
        J=J(2:end,:);
    end
end
J=reshape(J,[size(J,1),vecdim]);
end
