function [xdtvres,xdtvJ]=sys_cond_POEV1_selfphase_condition(point,dim)
%% obtain condition <x',v>=0
p1=setfield(point,'profile',point.profile(dim+(1:dim),:));  %#ok<*SFLD>
p2=setfield(point,'profile',point.profile(1:dim,:));
[xdtvres,pW]=dde_coll_profile_dot(p1,p2,'derivatives',[0,1]);
xdtvJ=p_axpy(0,point,[]);
xdtvJ.profile=cat(1,reshape(pW'*p1.profile(:),dim,[]),...
    reshape(pW*p2.profile(:),dim,[]));
end
