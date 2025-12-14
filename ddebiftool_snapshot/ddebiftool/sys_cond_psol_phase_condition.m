function [resph,Jph]=sys_cond_psol_phase_condition(point,pref,dim)
%% obtain condition <xref',x>=0
if nargin<3
    dim=size(point.profile,1);
end
pref.profile(dim+1:end,:)=0;
[resph,pW]=dde_coll_profile_dot(point,pref,'derivatives',[0,1]);
resph=resph+0.5*diff(sum(pref.profile(:,[end,1]).^2,1),[],2);
Jph=p_axpy(0,point,[]);
Jph.profile=reshape(pW*pref.profile(:),size(pref.profile));
end
