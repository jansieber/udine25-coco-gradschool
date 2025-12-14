function pt=dde_hcli_create(varargin)
%% create hcli point
%
% $Id: dde_hcli_create.m 308 2018-10-28 15:08:12Z jansieber $
%%
default={'kind','hcli','parameter',[],'mesh',[],'degree',[],'profile',[],'period',[],...
    'x1',[],'x2',[],'lambda_v',[],'lambda_w',[],'v',[],'w',[],'alpha',[],'epsilon',[],'nvec',[]};
pt=dde_set_options(default,varargin,'pass_on','point');
pt.kind='hcli';
if ~isempty(pt.profile)
    pt=dde_coll_check(pt);
end
% extra line to normalize alpha's
if norm(pt.alpha)~=0
    pt.alpha=pt.alpha/norm(pt.alpha);
end
end