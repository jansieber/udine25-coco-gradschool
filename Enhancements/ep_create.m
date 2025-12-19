function pt=ep_create(varargin)
%% create ep struct for bifurcation diagram outout
%
%%
default={'kind','ep','parameter',[],'x',[]};
[pt,dum,userdefined]=sco_set_options(default,varargin,'pass_on','point'); %#ok<ASGLU>
pt.kind='ep';
