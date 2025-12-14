function [freq,stst,ind]=dde_stst_critical_freq(funcs,stst,kind,stbmth,varargin)
default={'excludefreqs',[],'includehopf',false,'eigenvalue',[],...
    'evfield','l1','closest',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
ind=NaN;
%% for fold return 0 freq
if strcmp(kind,'fold')
    freq=0;
    return
end
%% if frequency prescribed return frequency
if ~isempty(options.eigenvalue)
    freq=abs(imag(options.eigenvalue));
    return
end
%% create full list of excluded frequencies
excludefreqs=abs(options.excludefreqs(:).');
if isfield(stst,'omega') && ~isempty(stst.omega) && ~options.includehopf ...
        && ~isnan(stst.omega)
    excludefreqs=[abs(stst.omega),excludefreqs];
end
%% if eigenvalues for point not yet computed, do so
if isempty(stst.stability)
    if isempty(funcs)
        error('dde_stst_critical_freq:arguments',...
            'dde_stst_critical _freq: eigenvectors not present in stst and r.h.s. not provided');
    end
    stst.stability=p_stabil(funcs,stst,stbmth,pass_on{:});
end
%% select eigenvalue field
evals=stst.stability.(options.evfield);
ind_im_ge0=1:length(evals);
%% remove frequencies to be excluded
if ~isempty(excludefreqs)
    excludefreqs=[excludefreqs,-excludefreqs(excludefreqs>0)];
    ind=dde_match_complex(1i*excludefreqs(:),evals(ind_im_ge0));
    ind_im_ge0(ind)=[];
end
ind_im_ge0=ind_im_ge0(imag(evals(ind_im_ge0))>=0);
if isempty(ind_im_ge0)
    error('dde_stst_critical_freq:eigenvalues',['dde_stst_critical_freq: ',...
        'no eigenvalues found.']);
end
%% select eigenvalue closest to imaginary axis or to options.closest
evals=evals(ind_im_ge0);
if isempty(options.closest) || ~isnumeric(options.closest)
    [i1,i2]=min(abs(real(evals))); %#ok<ASGLU>
    freq=abs(imag(evals(i2)));
else
    [i1,i2]=min(abs(evals-options.closest)); %#ok<ASGLU>
    freq=abs(imag(evals(i2)));
end
ind=ind_im_ge0(i2);
end
