function [funcs,branch,nremove]=PsolFromPsolbif(funcs,branch)
%% reduce branch of Torusbifurcations (period doublings) or POfolds to Psol
% $Id: PsolFromPsolbif.m 375 2019-09-14 14:26:50Z jansieber $
%%
nremove=[];
if isfield(funcs,'kind')
    switch funcs.kind
        case 'POfold'
            nremove=1;
        case 'PD'
            nremove=-1;
        case 'torus'
            omega=funcs.get_comp(branch,'omega');
            nremove=exp([-1,1]*1i*pi*omega);
    end
    branch.point=funcs.get_comp(branch.point,'solution');
    branch.parameter.free=intersect(branch.parameter.free,1:funcs.ip.nuserpar);
    branch.method=funcs.usermethod;
    funcs=funcs.userfuncs;
end
end