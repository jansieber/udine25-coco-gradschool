function [res,J]=rot_cond(point,funcs,pref)
%% extra phase condition needed to fix rotation phase
%
% $Id: rot_cond.m 357 2019-06-30 23:59:31Z jansieber $
%
A=funcs.rotation;
if isfield(point,'x')
    res=pref.x'*A*point.x;
    J=p_axpy(0,point,[]);
    J.x=A'*pref.x;
elseif isfield(point,'profile')
    pref.profile=A*pref.profile;
    pref.period=0;
    [res,J]=p_dot(point,pref,'free_par_ind',[]);
else
    error('rot_cond:type','rot_cond: type %s not supported',point.kind);
end
end
