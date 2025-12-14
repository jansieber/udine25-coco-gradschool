function [trfuncs,extra_freepar,initfuncs]=set_PDfuncs(funcs,point,method,varargin)
%% set up extended systems, numbering and values of additional parameters and artificial delays
%% process options
default={'biftype','PD'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
[trfuncs,extra_freepar,initfuncs]=set_torusfuncs(funcs,point,method,...
    'biftype',options.biftype,pass_on{:});
end