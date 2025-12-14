function [r,J,ieq]=dde_stst_nullphase_cond(pt,pref,free_par,method)
[r,J,ieq]=deal(zeros(0,1),repmat(pt,0,1),zeros(0,1));
if isempty(pref) || ...
        (nargin>3 && isfield(method,'phase_condition') && ~method.phase_condition)
    return
end
r=imag(pref.v'*pt.v);
J=p_axpy(0,pt,[]);
J.v=1i*pref.v;
ieq=1;
end
