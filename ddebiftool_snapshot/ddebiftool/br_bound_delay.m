function branch=br_bound_delay(funcs,branch)
%% bound delays from below by zero 
% if not yet done, delays are constant and requested by method.continuation
if funcs.tp_del || ...
        (isfield(branch.method.continuation,'permit_negative_delay')&&...
        ~branch.method.continuation.permit_negative_delay)
    return
end
tau=funcs.sys_tau();
[settau,ibd]=intersect(branch.parameter.min_bound(:,1),tau);
branch.parameter.min_bound(ibd,2)=max(branch.parameter.min_bound(ibd,2),0);
tau=setdiff(tau,settau);
nbd=size(branch.parameter.min_bound,1);
for j=1:length(tau)
    branch.parameter.min_bound(nbd+j,:)=[tau(j) 0];
end
end
