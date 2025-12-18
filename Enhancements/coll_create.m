function pt=coll_create(varargin)
%% create coll structure with guaranteed fields
%%
default={'kind','coll','parameter',[],'t',[],'degree',[],'x',[],...
    'T',[],'T0',0};
[pt,dum,userdefined]=sco_set_options(default,varargin,'pass_on','point'); %#ok<ASGLU>
pt.kind='coll';
end