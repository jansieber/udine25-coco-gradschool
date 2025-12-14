function [r,J]=sys_cond_RWFold(p,A,pref)
%% constraints used for extended DDE in fold continuation of relative equilibria
%
% fix phase and norm of derivative vector
J=p_axpy(0,p,[]);
r=pref.x'*A*p.v;
J.v=A'*pref.x;
end
