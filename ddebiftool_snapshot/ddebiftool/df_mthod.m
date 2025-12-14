function method=df_mthod(varargin)
%% create methed structure with default parameters for continuation, solving and stability
% function method=df_mthod(kind,flag_newhheur)
% INPUT:
%	kind the kind of default method wanted
%   opt: discretization (optional, default: 2) boolean: use Chebyshev approx
%                     Only used if kind==stst, alternatives: 0=bdf, 1=mxo
%                     (Verheyden)
% OUTPUT:
%	method default method
% COMMENT:
%       The order of the LMS method used in the computation of the
%       characteristic roots is hard coded in this file.
%  
% Update on 05/03/2007 ("flag_newhheur" <=> (imag(method.stability.lms_parameter_rho)~=0) )
%
% $Id: df_mthod.m 369 2019-08-27 00:07:02Z jansieber $
%
%% correction of argument list for compatibility
% (formerly first arg funcs no longer needed)
flag_newhheur=2; % compatibility to Verheyden format (but switching away from Verheyden's method)
if length(varargin)==1
    kind=varargin{1};
    pass_on=varargin(2:end);
elseif length(varargin)==2 && isstruct(varargin{1})
    kind=varargin{2};
    pass_on=varargin(3:end);
elseif length(varargin)==2 && ischar(varargin{1})
    kind=varargin{1};
    flag_newhheur=varargin{2};
    pass_on=varargin(3:end);
elseif length(varargin)>=3
    kind=varargin{2};
    flag_newhheur=varargin{3};
    pass_on=varargin(4:end);
else
    error('df_mthod:args','df_mthod: calling arguments unexpected');
end
try
    method=dde_apply({'dde_',kind,'_method'},pass_on{:});
    return
catch ME
    if ~strcmp(ME.identifier,'dde_apply:function')
        rethrow(ME);
    end
end
%% continuation parameters
method.continuation.steplength_condition=1; % use steplength condition
method.continuation.plot=1; % plot new points
method.continuation.prediction=1; % use secant prediction
method.continuation.steplength_growth_factor=1.2; % grow steplength with factor 1.2
method.continuation.steplength_minimum=0; % minimal stepsize
method.continuation.plot_progress=1; % plot progress gradually
method.continuation.plot_measure=[]; % use default plot measures
method.continuation.print_progress=0; % print out information during continuation
method.continuation.warnings=1; % print warnings (eg boundary reached,angles sharp)
method.continuation.warn_angle=0.5; % warning cos(angle) between correction and prediction
method.continuation.halt_before_reject=0; % rejection of points allowed
method.continuation.minimal_angle=-1; % minimal cos(angle) between prediction and correction
method.continuation.use_tangent=false; % use tangent as pseudoarclength condition
method.continuation.permit_negative_delay=false; % stop when delays become negative
method.continuation=br_add_stop(method.continuation); % functions indicating stopping criteria
%% Bifurcation detection
method.bifurcation.radial_tolerance_factor = 0.25;
method.bifurcation.minimal_real_part = -0.1;
method.bifurcation.correction_tolerance = 1e-6;
method.bifurcation.singular_tolerance = 1e-12;
method.bifurcation.secant_iterations = 30;
method.bifurcation.secant_tolerance = 1e-6;
method.bifurcation.print = 1;
method.bifurcation.imagthreshold = 1e-6;
method.bifurcation.monitor_eigenvalues = 0;
method.bifurcation.plot_testfunctions = 0;
method.bifurcation.pause_on_bifurcation= 0;
%% Newton iteration
method.point.newton_max_iterations=5;
if strcmp(kind,'hcli')
    method.point.newton_max_iterations=10;
    method.point.bctype='eigval'; % set 'subspace' to use subspace tracking (not yet implemented)
end
method.point.newton_nmon_iterations=1; % permit iterations to grow residual initially
method.point.preprocess='';            % prepare point before entering Newton ieration
method.point.postprocess='';           % postprocess point before entering Newton ieration (e.g. adjust phase etc)
method.point.extra_condition=true;     % add extra conditions (sys_cond's)
method.point.print_residual_info=0;    % print residual during Newton iteration
method.point.jacobian_nonsquare=false; % warn if Jacobian is non-square
method.point.remesh=false;             % perform remeshing controlled by other flags
method.point.extracolumns=[];          % add extra columns for overdetermined systems 
method.point.skip_jacobian=0;          % if new residual<skip_jacobian*old_residual, do not recompute Jacobian
switch kind
    case 'stst'
        % stst
        method.point.halting_accuracy=1e-10; % stop iteration when residual<halting
        method.point.minimal_accuracy=1e-8;  % consider iteration failed if final residual>minimal
    case 'fold'
        method.point.norm=true;              % add norm of eigenvector to equation
        method.point.halting_accuracy=1e-9;  % stop iteration when residual<halting
        method.point.minimal_accuracy=1e-7;  % consider iteration failed if final residual>minimal
    case 'hopf'
        method.point.norm=true;              % add norm of eigenvector to equation
        method.point.phase_condition=true;   % add phase condition for complex eigenvector
        method.point.preprocess='';%'dde_jac2square_preprocess';
        method.point.postprocess='';%'dde_hopf_postprocess';
        method.point.halting_accuracy=1e-9;  % stop iteration when residual<halting
        method.point.minimal_accuracy=1e-7;  % consider iteration failed if final residual>minimal
    case {'psol','hcli'}
        % point
        method.point.halting_accuracy=1e-8; % stop iteration when residual<halting
        method.point.minimal_accuracy=1e-6; % consider iteration failed if final residual>minimal
        method.point.phase_condition=1;     % add phase condition wrt reference for profile
        method.point.collocation_parameters=[]; % which collocation points in interval to use
        method.point.adapt_mesh_before_correct=0; % call p_remesh before correction every ? step
        method.point.adapt_mesh_after_correct=3; % call p_remesh after correction every ? step and correct again
        method.point.remesh=true;           % switch remeshing on/off
        method.point.matrix='full';
        method.point.maxsize=2^23;
        if strcmp(kind,'hcli')              % approximate integrals in invariant subspace projections
            [lgt,lgw]=legpts(12,[0,1]);     % using Gauss-Legendre formula
            method.point.lgpts=lgt;         % with lgt points and lgw weights
            method.point.lgw=lgw;
        end
    otherwise
        [dum1,dum2,kindparent]=feval(['dde_',kind,'_vars']); %#ok<ASGLU>
        if ~isempty(kindparent)
            method=df_mthod(kindparent,flag_newhheur,pass_on{:});
            return
        else
            warning('df_mthod: kind %s not recognized',kind);
        end
end
%% Stability/spectrum (cheb is default for stst, fold, hopf)
discretizations={'bdf','mxo','cheb'}; % methods available for stst, hopf, fold
method.stability.delay_accuracy=-1e-8;
method.stability.max_number_of_eigenvalues=20;
method.stability.root_accuracy=1e-6;
method.stability.fill=[]; % permit stability fields l0, l1 to have different field lengths (NaN fills with NaN)
switch kind
    case {'stst','hopf','fold'}
        if isnumeric(flag_newhheur)
            discretization=discretizations{flag_newhheur+1};
        elseif ischar(flag_newhheur)
            discretization=flag_newhheur;
        end
        method.stability.discretization=discretization;
        method.stability.delay_accuracy=-1e-8;
        method.stability.root_accuracy=1e-6;
        switch discretization
            case{'bdf','mxo'}
                order=4;
                delta_region=0.1;
                [alpha,beta]=dde_stst_time_lms(discretization,4);
                method.stability.lms_parameter_alpha=alpha;
                method.stability.lms_parameter_beta=beta;
                switch discretization
                    case 'mxo'
                        method.stability.lms_parameter_rho=dde_stst_mxo_ellipse(order,delta_region);
                        method.stability.newheuristics_tests=2000;
                    case 'bdf'
                        method.stability.lms_parameter_rho=dde_stst_bdf_safety(alpha,beta,0.01,0.01);
                end
                method.stability.interpolation_order=order;
                method.stability.minimal_time_step=0.01;
                method.stability.maximal_time_step=0.1;
                method.stability.max_newton_iterations=6;
                method.stability.remove_unconverged_roots=1;
                method.stability.minimal_real_part=[];
            case 'cheb'
                method.stability.max_number_of_eigenvalues=20;
                method.stability.min_number_of_eigenvalues=[];
                method.stability.ncheb=[]; % initial degree
                method.stability.inisize=[];
                method.stability.maxsize=200; % maximal size of matrix created
                method.stability.minimal_real_part=-1;       
                method.stability.scale_w=true;
                method.stability.nearest=[];
                method.stability.tauscal_limit=1e-5;
                method.stability.discard_accuracy_factor=1e5; % if finite, roots with err>root_accuracy*this will be discarded
        end
    case 'psol'   % multipliers
        method.stability.eigmatrix='full';
        method.stability.eigmaxsize=2^23; % maximum
        method.stability.closest=[];
        method.stability.fill=[];
        method.stability.minimal_modulus=0;
        method.stability.root_accuracy=1e-6;
        method.stability.collocation_parameters=[];
        method.stability.geteigenfuncs=false;
    case 'hcli'
        method=rmfield(method,'stability');
end

