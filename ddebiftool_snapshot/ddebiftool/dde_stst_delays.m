function [taus,xtau_values]=dde_stst_delays(funcs,pt,varargin)
%% compute all or requested delays in stst, including first delay=0
%%
default={'d_nr',[],'repeat',true,'max_xderiv',0};
options=dde_set_options(default,varargin,'pass_on');
xvec=cat(3,pt.x);
pvec=cat(3,pt.parameter);
ntaum1=dde_num_delays(funcs); % derivatives of arguments needed?
d=ntaum1+1; % number of delays incl 0
xtau_values=repmat(xvec,1,d,1,options.max_xderiv+1);
if options.max_xderiv>0
    xtau_values(:,:,:,2:end)=0;
    xtau_values=permute(xtau_values,[1,2,4,3]);
end
inc0=true;
taus=dde_taufunvec(funcs,xtau_values,pvec,options.repeat,inc0);
if ~isempty(options.d_nr)
    taus=taus(options.d_nr+1,:);
end
end
