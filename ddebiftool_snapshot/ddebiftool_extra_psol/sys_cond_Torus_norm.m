function [r,J]=sys_cond_Torus_norm(point,dim,varargin)
default={'res',1};
options=dde_set_options(default,varargin,'pass_on');
%% keep Floquet mode at unit length
Jtemplate=p_axpy(0,point,[]);
irgu=dim+(1:dim);
irgv=dim*2+(1:dim);
irguv=[irgu,irgv];
uvpoint=Jtemplate;
uvpoint.profile(irguv,:)=point.profile(irguv,:);
[utuvtv,W1]=dde_coll_profile_dot(uvpoint,uvpoint);
J=Jtemplate;
J.profile(:)=2*W1*uvpoint.profile(:);
r=utuvtv-options.res;
end
