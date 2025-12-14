function [points,suc]=dde_correct_resonance(trfuncs,trbranch,points,varargin)
%% correct located resonant points along torus bifurcation with p_correc
%% update branch method parameters with possible optional arguments
free_par_ind=trbranch.parameter.free;
free_par_ind(free_par_ind==trfuncs.ip.omega)=[];
trbranch=replace_branch_pars(trbranch,free_par_ind,varargin);
%% fix rotation number for each point and apply p_correc
mth=trbranch.method.point;
np=length(points);
suc=NaN(1,np);
for i=1:np
    [points(i),suc(i)]=p_correc(trfuncs,points(i),free_par_ind,[],...
        mth,1,points(i));
    if trbranch.method.continuation.print_progress>0
        fprintf('i=%d of %d, omega=%s, suc=%d\n',i,np,...
            rats(points(i).parameter(trfuncs.ip.omega)),suc(i));
    end
end
end
