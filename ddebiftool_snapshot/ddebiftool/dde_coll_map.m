function [coll,extmesh]=dde_coll_map(funcs,pt,varargin)
%% Map values on mesh in point structure to collocation points
%
%% Inputs:
% * funcs: problem definition  (sys_tau, sys_ntau, tpdel needed)
% * pt: point structure (eg, psol or hcli) with collocation mesh, degree and
% profile, and period
% * mesh evaluation function
%% Optional (name-value pairs):
% * wrapJ (true): wrap evaluation of points around periodically and adjust
% Jacobian correspondingly (for psol)
% * c_is_tvals: (false) whether collocation points are given explicitly in
% relation to the entire period
% * c (1 x pt.degree or 1 x neqs array): if c_is_tvals is false then tese
% are the points where collocation interpolation is evaluated. Otherwise,
% these are the collocation points in one subinterval (unscaled). If c is
% empty then the Gauss-Legendre points are used.
% * nderivs (default 1) number of derivatives to be computed
%
%% Outputs:
% * y n x (nderivs+1) x neqs array interpolation values
% * W (ntau+1) x (nderivs+1) cell array of (n*neqs) x numel(pt.profile)
% sparse Jacobian matrices
% * tc (1 x neqs) mesh of collocation points used (equal c if c is
% nonempty and c_is_tvals is true)
% * tau neqs x (ntau+1) array of delays
% * extmesh (if wrapJ is false, this contains the mesh points for all
%  all columns of y{k,d}.
%
% $Id: dde_coll_map.m 369 2019-08-27 00:07:02Z jansieber $
%%
default={'wrapJ',true,'max_xderiv',1,'scal_deriv',true};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[tc,submesh_limit]=coll_map_times(funcs,pt,pass_on{:});
[taus,y]=dde_coll_delays(funcs,pt,'t',tc,...
    'submesh_limit',submesh_limit,pass_on{:},'repeat',false,...
    'max_xderiv',options.max_xderiv);
%% generate W, W', W^(k-1) for each delay
% W{t_i}) is W for delay t_i-1 (t_i=1 corresponds to tau=0), W{t_i}{2}
% is W' at delay t_i-1.
taus=taus.';
d=size(taus,2);
W=cell(d,options.max_xderiv+1);
t_eval=tc(ones(d,1),:).'-taus/pt.period;
T_scal=pt.period^double(options.scal_deriv);
for k=0:size(y,3)-1
    [W(:,k+1),extmesh]=pt_eval(pt,t_eval,k,options.scal_deriv,'kron',true,...
        'wrapJ',options.wrapJ,'submesh_limit',submesh_limit,pass_on{:});
    y(:,:,k+1,:)=y(:,:,k+1,:)./T_scal^k;
end
max_rhs_xderiv=dde_num_delays(funcs,'max_rhs_xderiv');
coll=struct('mesh',tc,'profile',y,'jac',{W},'tau',taus/pt.period,...
    'max_rhs_xderiv',max_rhs_xderiv,'period',pt.period);
end
%%
function [W,extmesh]=pt_eval(pt,x,order,scal_deriv,varargin)
n=size(pt.profile,1);
[neqs,ntau_i]=size(x);
[Wc,dum,ptext]=dde_coll_wrap_eva(pt,x(:).','diff',order,varargin{:},...
    'output','matrix'); %#ok<ASGLU>
extmesh=ptext.mesh;
scal=pt.period.^(double(scal_deriv).*order);
W=mat2cell(Wc/scal,n*neqs*ones(ntau_i,1),numel(ptext.profile));
end

function [tc,submesh_limit] = coll_map_times(funcs,pt,varargin)
%% Determine where residuals are computed: mesh tc 
% only compute res & jac at specified points?
default={'c_is_tvals',false,{'c','collocation_parameters'},[],'submesh_limit',0};
options=dde_set_options(default,varargin,'pass_on');
submesh_limit=options.submesh_limit;
if ~options.c_is_tvals
    if isempty(options.c) || ischar(options.c)
        options.c=dde_coll_set_grid('collocation',pt.degree,'type',options.c,...
            'lhs_num',funcs.lhs_matrixfun(size(pt.profile,1)));
    end
    if 1==max(options.c)
        submesh_limit=1;
    end
    tc=dde_coll_meshfill(pt.mesh(1:pt.degree:end),pt.degree,...
        'purpose','collocation','grid',options.c);
else
    tc=options.c(:)';
end
end