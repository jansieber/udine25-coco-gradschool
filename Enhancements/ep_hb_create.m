function pt=ep_hb_create(varargin)
%% create Hopf point 
default={'kind','hopf','parameter',[],'x',[],'v',[],'k'};
[pt,dum,userdefined]=dde_set_options(default,varargin,'pass_on','point'); %#ok<ASGLU>
pt.kind='hopf';
if ~isempty(pt.x) && isempty(pt.v)
    pt.v=NaN(size(pt.x));
end
if ~isempty(pt.x) && isempty(pt.k)
    pt.k=NaN;
end
end
