%% Create initial q-periodic orbit near psol with Floquet multiplier exp(1i*2*pi*p/q)
%%
function [perq,tangent]=dde_psol_from_psol(point,varargin)
%% Inputs
% 
% * |point|: |'psol'| periodic orbits from which one wants to branch off
% 
% Important optional inputs (name-value pairs)
%
% * |funcs| (mandatory): structure with functions provided by user
% * |'radius'|: initial deviation along period-doubling eigenvector
% * |'resonance'|: (1x2 integers ,[p,q]) q-periodic orbit is constructed,
% default [1,2]
%
% All other name-value pairs are passed on to output branch.
%% Outputs
%
% * |perq|: periodic orbit with initial small period 2 deviation
% * |tangent|: approximate tangent
%
% $Id: dde_psol_from_psol.m 308 2018-10-28 15:08:12Z jansieber $
%
%%
default={'radius',0.01,'method',getfield(df_mthod('psol'),'stability'),'funcs',[],'resonance',[1,2],...
    'angle_tol',0.25,'rotation',[0,1],'initmethod','eig'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
% create branch perq of periodic solutions starting from an
% approximate period doubling num on a branch br of periodic orbits
% obtain eigenvector closest to exp(2*pi*1i*p/q
p_by_q=options.resonance(1)/options.resonance(2);
floquet=exp(2*1i*pi*p_by_q);
tangent=repmat(point,0,1);
%% find eigenvector for eigenvalue closest to floquet
if strcmp(options.initmethod,'eig')
    if  ~isfield(point,'stability') || ~isfield(point.stability,'eigenfuncs')...
            || isempty(point.stability.eigenfuncs)
        if isempty(options.funcs)
            warning('dde_psol_from_psol:arguments',...
                'dde_psol_from_psol: eigenvectors and equations not provided');
            perq=point;
            return
        end
        stability=p_stabil(options.funcs,point,options.method,'closest',floquet,...
            'geteigenfuncs',true,'max_number_of_eigenvalues',1,pass_on{:});
    else
        stability=point.stability;
    end
    [~,icrit]=min(abs(stability.mu-floquet));
    eigval=stability.mu(icrit);
    omega=atan2(imag(eigval),real(eigval))/pi/2;
    if abs(omega-p_by_q)>=options.angle_tol
        warning('dde_psol_from_psol:eigenvalue',...
            'dde_psol_from_psol: angle for eigenvalue closest to exp(2pi i %d/%d), %g+1i%g, differs by %g degree',...
            options.resonance(1),options.resonance(2),real(eigval),imag(eigval),(omega-p_by_q)*180);
        perq=point;
        return
    end
    cdeviation=stability.eigenfuncs(icrit).profile;
    basepoint=point;
else % initmethod 'svd'
    assert(~isempty(options.funcs)&&isfield(options.funcs,'ip'))
    ip=options.funcs.ip;
    cdeviation=point.profile(ip.dim+(1:ip.extdim),:);
    basepoint=options.funcs.get_comp(point,'solution');
end
%% use real part of eigenvector to construct small deviation from point
% call it perq. tangent equals eigenvector
q=options.resonance(2);
crot=exp(1i*2*pi*options.rotation(1)/options.rotation(2));
rotdeviation=crot*cdeviation;
[nx,nt]=size(rotdeviation);
cdevrep=repmat(rotdeviation,1,1,q).*repmat(reshape(floquet.^(0:q-1),1,1,q),nx,nt,1);
cdeviationlong=[rotdeviation(:,1),reshape(cdevrep(:,2:end,:),nx,(nt-1)*q)];
deviation=real(cdeviationlong);
deviation=deviation/max(abs(deviation(:)));
tmeshrep=(repmat((0:q-1),nt,1)+repmat(basepoint.mesh',1,q))/q;
tmesh=[tmeshrep(1),reshape(tmeshrep(2:end,:),1,[])];
perq=dde_psol_create('point',basepoint,...
    'mesh',tmesh,...
    'period',q*basepoint.period,...
    'profile',[basepoint.profile,repmat(basepoint.profile(:,2:end),1,q-1)]+options.radius*deviation);
tangent=dde_psol_create('point',perq,...
    'parameter',0*perq.parameter,...
    'period',0,...
    'profile',deviation);
end
