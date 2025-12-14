%% normalization of nullvector
function [r,J,ieq]=dde_stst_nullnorm_cond(pt,free_par,method,ref)
[r,J,ieq]=deal(zeros(0,1),repmat(pt,0,1),zeros(0,1));
if nargin>2 && isfield(method,'norm') && ~method.norm
    return
end
if nargin<4
    ref=1;
end
r=pt.v'*pt.v+[0,pt.parameter(free_par)]*[0,pt.parameter(free_par)]'-ref;
J=p_axpy(0,pt,[]);
J.v=2*pt.v;
J.parameter(free_par)=2*pt.parameter(free_par);
ieq=1;
end