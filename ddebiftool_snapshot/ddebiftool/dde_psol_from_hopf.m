function [psol,tangent]=dde_psol_from_hopf(point,varargin)
%% Create small amplitude harmonic oscillation around Hopf point
%
% Output is initial periodic orbit and approximate tangent
% $Id: dde_psol_from_hopf.m 369 2019-08-27 00:07:02Z jansieber $
%%
default={'radius',1e-3,'degree',3,'intervals',20,...
    'submesh',[],'collocation_parameters',[],'funcs',[],...
    'harmonic_from_v',@(v,t)imag(v*exp(2*pi*1i*t))};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
if ~strcmp(point.kind,'hopf')
    point=dde_apply({'dde_hopf_from_',point.kind,''},point,pass_on{:});
end
if point.omega==0
    warning('dde_psol_from_hopf:zerofrequency',...
        'zero frequency in Hopf point, period set to Inf');
end
if length(options.intervals)==1
    coarsemesh=linspace(0,1,options.intervals+1);
else
    coarsemesh=options.intervals;
end
force_smooth=ischar(options.collocation_parameters) &&...
    strcmp(options.collocation_parameters,'force_smooth');
if force_smooth
    options.submesh='force_smooth';
end
pmesh=dde_coll_meshfill(coarsemesh,options.degree,'grid',options.submesh,...
    'lhs_num',loc_lhs_test(options.funcs,point.x));
psol=dde_psol_create('point',point,...
    'mesh',pmesh,...
    'degree',options.degree,...
    'period',2*pi/abs(point.omega),...
    'profile',repmat(point.x,1,length(pmesh)));
%% harmonic profile
harmonic=options.harmonic_from_v(point.v,psol.mesh);
%% initial solution has small amplitude harmonic aroun equilibrium
psol.profile=psol.profile+options.radius*harmonic;
tangent=dde_psol_create('point',psol,...
    'period',0,...
    'profile',harmonic,...
    'parameter',0*psol.parameter);
end
function mat=loc_lhs_test(funcs,x)
if isempty(funcs)
    mat=1;
else
    mat=funcs.lhs_matrixfun(size(x,1));
end
end