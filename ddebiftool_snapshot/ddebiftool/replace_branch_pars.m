function branch=replace_branch_pars(branch,contpar,pass_on)
%% pass on optional arguments to branch structure
%
% branch=replace_branch_pars(branch,contpar,pass_on)
% 
% input
% branch: given branch structure to be amended
% contpar: continuation parameter, prepended to free parameters already
%          present in branch (if length==1) or replacing
%          branch.parameter.free
% pass_on: cell array containing name-value pairs for fields of
%          branch.method.continuation, branch.method.point, branch.method.stability and
%          branch.parameter to be replaced
%          may also include name-value pair 'branchtype',kind, where kind
%          is the new type of branch. In this case defaults for the new
%          kind are set before the old values are copied over
%
% output: branch with amended fields
%
% $Id: replace_branch_pars.m 339 2019-05-09 19:47:01Z jansieber $
%
%%
default={'branchtype',[]};
[options,pass_on]=dde_set_options(default,pass_on,'pass_on');
branch.method=replace_method_pars(branch.method,options.branchtype,pass_on{:});
if isempty(contpar)
    % if no continuation parameters are given use free parameters of branch
    branch.parameter.free=branch.parameter.free;
else
    branch.parameter.free=contpar(:)';
end
[sel,other]=dde_options_filter('parameter',pass_on);
branch.parameter=dde_set_options(branch.parameter,sel,'pass_on');
branch.parameter=dde_set_options(branch.parameter,other,'pass_on');
end
