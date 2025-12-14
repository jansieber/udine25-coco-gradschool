function dfun=nmfm_dev_dvlam(fun,ord,dv,dlam)
%% Directional deriv of history function wrt v and lam of order ord
% in direction dv (size(v)), dlam (size(lambda))
%%
if ord==0
    dfun=fun;
    return
end
n=size(fun.v,1);
dlamkm1=dlam.^(ord-1);
dlamk=dlam.^ord;
v1=ord*dv.*dlamkm1(ones(n,1),:);
t1=fun.t+ord-1;
t2=fun.t+ord;
v2=fun.v.*dlamk(ones(n,1),:);
dfun=nmfm_dev_fun([v1,v2],'lambda',[fun.lambda,fun.lambda],'t',[t1,t2]);
end
