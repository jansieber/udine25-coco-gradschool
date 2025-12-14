function hopf=dde_hopf_from_stst(stst,varargin)
%% convert steady state stst into approximate hopf point (not corrected)
%
% optional arguments: 
%
% * 'excludefreqs': list of frequencies to exclude
% * 'funcs': r.h.s. functions structure to compute eigenvectors if needed
% * 'method': parameters for stability computation
% $Id: dde_hopf_from_stst.m 315 2019-01-29 19:42:21Z jansieber $
%%
default={'funcs',[],'stability',getfield(df_mthod('stst'),'stability'),...
    'method',getfield(df_mthod('hopf'),'point'),...
    'recompute',false,'v_scal',@(p,pref)real(p.v'*p.v)};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
hopf=dde_hopf_create('point',stst,'stability',stst.stability);
if ~options.recompute && isfield(stst,'nvec') && isfield(stst,'flag') && ...
        ~isempty(stst.nvec) && strcmp(stst.flag,'hopf')
    %% Hopf eigenvector and value has been computed during normal form computations
    hopf.nvec=stst.nvec;
    hopf.omega=abs(hopf.nvec.omega);
    hopf.v=hopf.nvec.q;
    hopf.v=hopf.v/norm(hopf.v);
    hopf.nmfm=stst.nmfm;
    return
end
%% find corresponding eigenvector
if isstruct(hopf.stability) && isfield(hopf.stability,'v') && options.recompute
    hopf.stability=[];
end
if (isempty(hopf.stability) || ~isfield(hopf.stability,'v')) && isempty(options.funcs)
    error('dde_hopf_from_stst:arguments',...
        'dde_hopf_from_stst: eigenvectors not present in stst and r.h.s. not provided');
end
[freq,hopf,ind]=dde_stst_critical_freq(options.funcs,hopf,'hopf',options.stability,pass_on{:});
if ~isempty(options.funcs)
    hopf=dde_stst_critical_space(options.funcs,hopf,'hopf',freq,options.method,'output','pini','kind','hopf',pass_on{:});
else
    hopf.omega=freq;
    hopf.v=hopf.stability.v(:,ind);
end
try
    rs=options.v_scal(hopf);
catch
    rs=options.v_scal(hopf,hopf);
end
hopf.v=hopf.v/sqrt(rs);
end