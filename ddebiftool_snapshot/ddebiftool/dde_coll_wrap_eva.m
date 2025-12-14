function [y,J,ptext]=dde_coll_wrap_eva(pt,x,varargin)
%% wrapper around dde_coll_eva permitting times outside mesh interval
% treatment depends on option wrapJ: if wrapJ then times are taken modulo
% mesh interval (but including upper boundary). Otherwise, profile and mesh
% of point pt are extended periodically (up to boundary difference).
% Additional output is point structure with extended mesh and profile (for
% wrapJ==false). For wrapJ==true ptext equals pt.
default={'wrapJ',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
ptext=pt;
switch options.wrapJ
    case true
        x=dde_psol_wrap('mesh',pt.mesh,'t',x);
    case false
        [ptext.mesh,ptext.profile]=dde_apply({'dde_',pt.kind,'_extend_profile'},...
            pt,x,'profile',pt.profile,pass_on{:});
end
yJ=cell(1,2);
[yJ{:}]=dde_coll_eva(ptext,x,pass_on{:});
[y,J]=deal(yJ{:});
end
