function pt=ep_sn_create(varargin)
%% create fold point with empty stability and normal form information
%%
default={'kind','fold','parameter',[],'x',[],'v',[]};
[pt,dum,userdefined]=sco_set_options(default,varargin,'pass_on','point'); %#ok<ASGLU>
pt.kind='ep.SN';
if ~isempty(pt.x) && isempty(pt.v)
    pt.v=NaN(size(pt.x));
end
end
