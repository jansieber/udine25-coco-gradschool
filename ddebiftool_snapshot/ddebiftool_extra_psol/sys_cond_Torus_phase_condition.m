function [r,J]=sys_cond_Torus_phase_condition(point,pref,dim)
%% ensure that real and imaginary part of Floquet mode are orthogonal
% to reference
inner_matrix=kron([0,0,0;0,0,-1;0,1,0],speye(dim));
[r,Worth]=dde_coll_profile_dot(point,pref,'inner_matrix',inner_matrix);
J=p_axpy(0,point,[]);
J.profile(:)=Worth*pref.profile(:);
end