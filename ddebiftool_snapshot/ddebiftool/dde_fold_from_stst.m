function fold=dde_fold_from_stst(stst,varargin)
%% convert steady state stst into approximate fold point (not corrected)
%
% optional arguments: 
%
% * 'funcs': r.h.s. functions structure to compute char. matrix if needed
%
% $Id: dde_fold_from_stst.m 308 2018-10-28 15:08:12Z jansieber $
%%
default={'funcs',[],'recompute',false,'v_scal',@(p,pref)real(p.v'*p.v),...
    'free_par',zeros(0,1),'method',getfield(df_mthod('fold'),'point')};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
fold=dde_fold_create('point',stst,'stability',stst.stability);
if  ~options.recompute  && isfield(stst,'nvec') && isfield(stst,'flag') && ...
        ~isempty(stst.nvec) && strcmp(stst.flag,'fold')
    %% Fold eigenvector has been computed during normal form computations
    fold.nvec=stst.nvec;
    fold.v=fold.nvec.q;
    fold.v=fold.v/norm(fold.v);
    fold.nmfm=stst.nmfm;
    return
end
if isstruct(fold.stability) && isfield(fold.stability,'v') && options.recompute
    fold.stability=[];
end
%% check if eigenvectors are provided, 
% if not, compute minimal eigenvalue of charactaristic matrix
if (isempty(fold.stability) || ~isfield(fold.stability,'v')) && isempty(options.funcs)
    error('dde_fold_from_stst:arguments',...
        'dde_fold_from_stst: eigenvectors not present in stst and r.h.s. not provided');
end
if isempty(fold.stability) || ~isfield(fold.stability,'v') || ...
        (isfield(options.funcs,'ip') && isfield(options.funcs.ip,'nullparind'))
    fold=dde_stst_critical_space(options.funcs,fold,'fold',0,options.method,'output','pini',...
        'kind','fold','free_par',options.free_par,pass_on{:});
else
   [i1,i2]=min(abs(fold.stability.l1)); %#ok<ASGLU> 
   fold.v=real(fold.stability.v(:,i2));
end
try
    rs=options.v_scal(fold);
catch
    rs=options.v_scal(fold,fold);
end
fold.v=fold.v/sqrt(rs);
end

