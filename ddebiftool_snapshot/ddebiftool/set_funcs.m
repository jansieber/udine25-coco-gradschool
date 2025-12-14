function funcs=set_funcs(varargin)
%% fill in funcs structure with user-defined functions for use with DDE-Biftool functions
%
% possible named arguments: 'sys_rhs' (mandatory),'sys_ntau', 'sys_cond',
%  'sys_dirderi', sys_dirdtau, 'x_vectorized' (logical), 'p_vectorized'
%  (logical),'lhs_matrix' (l.h.s.matrix), 'sys_tau_seq' (hint if multiple
%  delays can be evaluated simultaneously), 'delayed_derivs' (up to which
%  order of derivatives sys_rhs, sys_dirderi depend on, adding an
%  additional dimension to xx if nonzero). No longer recommended:
%  'sys_deri', 'sys_dtau',
%
% Simplest Example
% for DDE x'=-x(t-tau)+a*x^3-x^5 with parameters [tau,a]:
% funcs=set_funcs(...
%   'sys_rhs',@(x,p)-x(1,2)+p(2)*x(1,1)^3-x(1,1)^5,...
%   'sys_tau',@()1)
% uses default mult_deriv for directional derivatives, has no
% vectorization.
% funcs=set_funcs(...
%   'sys_rhs',@(x,p)-x(1,2,:)+p(1,2,:).*x(1,1,:).^3-x(1,1,:).^5,...
%   'sys_tau',@()1,'x_vectorized', true, 'p_vectorized',true)
% is the same but has vectorization (will be faster, especially for
% periodic orbits)
%% Process options
defaults={...
    'sys_rhs',[],...             % user-defined r.h.s.
    'sys_ntau',@()0,...          % number of delays (SD-DDEs only)
    'sys_tau',@()[],...          % index of delays as parameters or dependence of delays
    'sys_cond',[],...            % extra conditions
    'sys_deri',[],...            % Jacobians and second-order derivatives
    'sys_dtau',[],...            % Jacobians and second-order derivatives for delay (SD-DDEs only)
    'sys_mfderi',{},...          % higher-order derivatives (up to 5) for constant delay stst normal forms
    'sys_dirderi',{},...         % directional derivatives of sys_rhs
    'sys_dirdtau',{},...         % directional derivatives of sys_tau
    'wrap_rhs',[],...            % r.h.s. wrapped into fun_vec for pseudo vectorization (if provided sys_ ... is ignoried)
    'wrap_dirderi',{},...        % r.h.s. derivatives wrapped into fun_vec for pseudo vectorization (if provided sys_ ... is ignoried)
    'wrap_tau',[],...            % delays wrapped into fun_vec for pseudo vectorization (if provided sys_ ... is ignoried)
    'wrap_dirdtau',{},...        % delay derivatives wrapped into fun_vec for pseudo vectorization (if provided sys_ ... is ignoried)
    'drhs_mf', [],...            % multi-directional and expandable derivative (constructed from wrap_dirderi, but overridable)
    'drhs_dir',[],...            % wrapped directional derivative arbitrary order (w option 'hjac',constructed from wrap_dirderi, but option overrides)
    'dtau_mf',[],...             % multi-directional and expandable derivative (constructed from wrap_dirdtau, but overridable)
    'dtau_dir',[],...            % wrapped directional derivative arbitrary order (w option 'hjac',constructed from wrap_dirdtau, but options overrides)
    'x_vectorized',false,...     % can r.h.s, tau and their Jacobians be called with x=n x (nd+1) x nvec?
    'p_vectorized',false,...     % --"-- with x=n x (nd+1) x nvec and p= 1 x np x nvec?
    'hjac',@(ord)eps^(1/(2+ord)),... % deviation to be used for finite differences in numerical Jacobians
    'sys_cond_reference',false,...% does sys_cond need second (reference) input?
    'lhs_matrix',NaN,...         % for DDAEs an optional lhs_matrix can be provided
    'delayed_derivs',0,...       % derivatives for arguments appended up to order (default 0)
    'dirderi_num',[],...         % hint if derivatives are provided (if only drhs_.. and dtau_.. are given, otherwise, determined from wrap_dirderi, wrap_dirdtau)
    'sys_tau_seq',[],...         % hint, replaced by {1,...,sys_ntau()} if empty: can sys_tau be called 
     ...                         % for several tau simultaneously?
     'xpattern',[]};             % optional hint: which delayed & non-delayed x's do rhs and tau depend on?
if isstruct(varargin{1})&& isfield(varargin{1},'sys_rhs')
    funcs=loc_copy_funcs(varargin{1},varargin{2:end});
    return
else
    funcs=dde_set_options(defaults,varargin);
end
sys_rhs_cell=false;
if isempty(funcs.sys_rhs) && isempty(funcs.wrap_rhs)
    file=exist('sys_rhs','file');
    if file==2
        funcs.sys_rhs=@sys_rhs;
    %else
    %    error('sys_rhs undefined');
    end
elseif iscell(funcs.sys_rhs)
    sys_rhs_cell=true;
    sys_rhs_fmt=funcs.sys_rhs(1:2);
    funcs.sys_rhs=@(x,p)dde2sys_rhs(0,funcs.sys_rhs{:},{},x,p);
end
if isempty(funcs.sys_cond)
    funcs.sys_cond=repmat(dde_sys_cond_create(),0,1);
elseif ~isstruct(funcs.sys_cond)
    funcs.sys_cond=dde_sys_cond_create('name','sys_cond','fun',funcs.sys_cond,...
        'reference',funcs.sys_cond_reference);
else
    funcs.sys_cond=reshape(funcs.sys_cond,1,[]);
end
%% convert user-provided (or NaN) lhs matrix into function
% we do not know the dimension of the problem if the user does not specify
% lhs_matrix. Hence, funcs.lhs_matrix is a function that takes as its
% argument the dimension of x.
if isnan(funcs.lhs_matrix)
    funcs.is_lhs_matrix_set=false;
else  % in this case the argument of lhs_matrix_fun is irrelevant
    funcs.is_lhs_matrix_set=true;
end
funcs.lhs_matrixfun=@(varargin)lhs_matrix(funcs.lhs_matrix,varargin{:});
funcs=rmfield(funcs,'lhs_matrix');
%% test for state-dependent delay
funcs.tp_del=true;
if ~iscell(funcs.sys_tau)
    try
        if ~isempty(funcs.wrap_tau)
            dummytau=funcs.wrap_tau();
        else
            dummytau=funcs.sys_tau();
        end
        funcs.tp_del=false;
        funcs.sys_ntau=@()length(dummytau);
    catch %#ok<CTCH>
        funcs.tp_del=true;
    end
end
%% pseudo-vectorize and make sure that the output has expected shape
isvec=[funcs.x_vectorized,funcs.p_vectorized];
[xfdims,xtaudims,pdims]=deal(2+double(funcs.delayed_derivs>0),2,2);
if isempty(funcs.wrap_rhs)
    funcs.wrap_rhs = @(x,p)fun_vec(funcs.sys_rhs,{1,[xfdims,pdims]},isvec,size(funcs.lhs_matrixfun(size(x,1)),1),x,p);
end
%% sys_dirderi not provided but sys_deri
if ~isempty(funcs.sys_deri) 
    funcs.wrap_deri_old= @(x,p,nx,np,v)fun_vec(...
    @(xa,pa)funcs.sys_deri(xa,pa,nx,np,[]),{2,[xfdims,pdims]},isvec([1,2]),[],x,p);
    if isempty(funcs.sys_dirderi) && isempty(funcs.wrap_dirderi)
        funcs.sys_dirderi{1}=...
            @(x,p,dx,dp)dde_dirderi_from_deri(x,p,dx,dp,funcs);
    end
end
%% rhs derivatives
if sys_rhs_cell && ~isempty(funcs.sys_dirderi)
    for i=1:length(funcs.sys_dirderi)
        funcs.sys_dirderi{i}=@(x,p,dx,dp)dde2sys_rhs(i,sys_rhs_fmt{:},funcs.sys_dirderi{i},{},x,p,dx,dp);
    end
end
if isempty(funcs.wrap_dirderi) && ~isempty(funcs.sys_dirderi)
    vecrhs=@(f){@(x,p,dx,dp)fun_vec(f,[xfdims,pdims,xfdims,pdims],[isvec,isvec],size(funcs.lhs_matrixfun(size(x,1)),1),x,p,dx,dp)};
    funcs.wrap_dirderi = cellfun(vecrhs,funcs.sys_dirderi);
end
%% create mixed higher-order, full and single-directional derivatives of rhs
% of arbitrary order
funcs=loc_funcs_add_mfderiv(funcs,'rhs');
%% funcs.sys_deri not provided (but possibly dirderi)
% combine Jacobians from directional derivatives
if isempty(funcs.sys_deri)
    funcs.wrap_deri_old=...
        @(x,p,nx,np,v)dde_gen_deriv(funcs.drhs_mf,x,p,nx,np,v);
end
%% calling sequence for delays
if ~funcs.tp_del
    taupar=funcs.sys_tau();
    funcs.wrap_tau = @(itau,x,p)loc_taufunc(0,itau,x,p,[],taupar);
    funcs.wrap_dirdtau = arrayfun(@(i){...
        @(itau,x,p,dx,dp)loc_taufunc(i,itau,x,p,dp,taupar)},...
        1:5);
    funcs.sys_tau_seq={1:funcs.sys_ntau()};
    funcs=loc_funcs_add_mfderiv(funcs,'tau');
    return
end
%% tau functions and dervatives of tau only for sd-DDEs
sys_tau_cell=false;
if iscell(funcs.sys_tau)
    sys_tau_cell=true;
    funcs.sys_tau=@(it,varargin)funcs.sys_tau{it}(varargin{:});
end
if sys_rhs_cell
    funcs.sys_tau=@(it,x,p)dde2sys_rhs(0,sys_rhs_fmt{:},funcs.sys_tau,{it},x,p);
end
if isempty(funcs.wrap_tau)
    funcs.wrap_tau = @(itau,x,p)fun_vec(@(xa,pa)funcs.sys_tau(itau,xa,pa),...
        [xtaudims,pdims],isvec,length(itau),x,p);
end
if isempty(funcs.sys_tau_seq)
    funcs.sys_tau_seq=num2cell(1:funcs.sys_ntau());
end
%% if sys_dtau is provided and vectorization is on, wrap output to ensure right format
%% sys_dirdtau not provided but sys_dtau
if ~isempty(funcs.sys_dtau)
    funcs.wrap_dtau_old=@(itau,x,p,nx,np)fun_vec(@(xa,pa)funcs.sys_dtau(xa,pa,nx,np),...
        [xtaudims,pdims],isvec,[],x,p);
    if isempty(funcs.sys_dirdtau) && isempty(funcs.wrap_dirdtau)
        funcs.sys_dirdtau{1}=...
            @(itau,x,p,dx,dp)dde_dirderi_from_deri(1,x,p,dx,dp,funcs,itau);
    end
end
%% delay derivatives
if sys_tau_cell
    for i=1:length(funcs.sys_dirdtau)
        funcs.sys_dirdtau{i}=@(it,varargin)funcs.sys_dirdtau{i}{it}(varargin{:});
    end
end
if sys_rhs_cell && ~isempty(funcs.sys_dirdtau)
    for i=1:length(funcs.sys_dirdtau)
        funcs.sys_dirdtau{i}=@(it,x,p,dx,dp)dde2sys_rhs(i,sys_rhs_fmt{:},funcs.sys_dirdtau{i},{it},x,p,dx,dp);
    end
end
if isempty(funcs.wrap_dirdtau) && ~isempty(funcs.sys_dirdtau)
    vectau=@(f){@(it,x,p,dx,dp)fun_vec(@(xa,pa,dxa,dpa)f(it,xa,pa,dxa,dpa),...
        [xtaudims,pdims,xtaudims,pdims],[isvec,isvec],size(funcs.lhs_matrixfun(size(x,1)),1),x,p,dx,dp)};
    funcs.wrap_dirdtau =cellfun(vectau,funcs.sys_dirdtau);
end
funcs=loc_funcs_add_mfderiv(funcs,'tau');
%% sys_dtau not provided (but possibly dirdtau)
if isempty(funcs.sys_dtau)
    funcs.wrap_dtau_old=@(itau,x,p,nx,np)dde_gen_deriv(...
        @(varargin)funcs.dtau_mf(itau,varargin),x,p,nx,np,[]);
end
end
%% uniform format for lhs_matrix
function L=lhs_matrix(arg,sz)
if isnan(arg)
    L=eye(sz(1));
elseif nargin>1
    L=arg;
else
    error('lhs_matrix:xunknown','lhs_matrix not given and no input provided');
end
end
%% append wrapper for derivatives and expansion to funcs
function out=loc_funcs_add_mfderiv(funcs,type)
out=funcs;
[directional,inc0]=deal(true);
switch type
    case 'rhs'
        fcn={funcs.wrap_rhs,funcs.wrap_dirderi};
        out=loc_check_add(out,'drhs_mf',...
            @(varargin)dde_mult_deriv(funcs,'rhs',fcn,~directional,varargin{:}));
        out=loc_check_add(out,'drhs_dir',...
            @(varargin)dde_mult_deriv(funcs,'rhs',fcn, directional,varargin{:}));
    case 'tau'
        fcn={@(x,p)dde_taufunvec(funcs,x,p,false,inc0),...
            arrayfun(@(ord){@(xa,pa,dx,dp)dde_all_dirdtau(funcs,inc0,ord,xa,pa,dx,dp)},...
            1:length(funcs.wrap_dirdtau))};
        out=loc_check_add(out,'dtau_mf',...
            @(varargin)dde_mult_deriv(funcs,'tau',fcn,~directional,varargin{:}));
        out=loc_check_add(out,'dtau_dir',...
            @(varargin)dde_mult_deriv(funcs,'tau',fcn, directional,varargin{:}));
end
out.dirderi_provided=@()loc_dirderi_provided(out);
end
%%
function funcs=loc_check_add(funcs,fname,value)
if ~isfield(funcs,fname) || isempty(funcs.(fname))
   funcs.(fname)=value;
end
end
%%
function tau=loc_taufunc(order,it,xx,p,dp,itau) %#ok<INUSL>
%% when converting DDE with parameter delays to sd-DDE, 
% this is the sys_tau and sys_dirdtau
%%
if order==0
    tau=p(1,itau(it),:);
elseif order==1
    tau=dp(1,itau(it),:);
else
    tau=0*dp(1,ones(1,length(it)),:);
end
end
%% 
function order=loc_dirderi_provided(funcs)
if ~isempty(funcs.dirderi_num)
    order=funcs.dirderi_num;
    return
end
order=length(funcs.wrap_dirderi);
if ~funcs.tp_del
    return
end
order=min(order,length(funcs.wrap_dirdtau));
end
function outfuncs=loc_copy_funcs(funcs,varargin)
default={'wrap_rhs',funcs.wrap_rhs,'wrap_dirderi',funcs.wrap_dirderi,...
    'wrap_tau',funcs.wrap_tau,'wrap_dirdtau',funcs.wrap_dirdtau,...
    'sys_ntau',funcs.sys_ntau,'sys_tau_seq',funcs.sys_tau_seq,...
    'sys_cond',funcs.sys_cond,'sys_cond_reference',funcs.sys_cond_reference,...
    'xpattern',funcs.xpattern,'hjac',funcs.hjac,...
    'delayed_derivs',funcs.delayed_derivs,...
    'x_vectorized',true,'p_vectorized',true};
outfuncs=set_funcs(default{:},varargin{:});
if funcs.is_lhs_matrix_set && ~outfuncs.is_lhs_matrix_set
    outfuncs.lhs_matrix_fun=funcs.lhs_matrix_fun;
end
end