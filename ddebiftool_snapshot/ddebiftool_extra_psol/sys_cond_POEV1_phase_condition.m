function [xdtvres,xdtvJ]=sys_cond_POEV1_phase_condition(point,pref,dim)
%% obtain condition <xref',v>=0
pref.profile(dim+(1:dim),:)=pref.profile(1:dim,:);
pref.profile(1:dim,:)=0;
[xdtvres,pW]=dde_coll_profile_dot(pref,point,'derivatives',[1,0]);
xdtvJ=p_axpy(0,point,[]);
xdtvJ.profile=reshape(pref.profile(:)'*pW,size(pref.profile));
end
