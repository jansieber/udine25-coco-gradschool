function method=replace_method_pars(inp,branchtype,varargin)
%% overwrite method parameters with argument list
%
%%
if ~isempty(branchtype)
    try
        method=df_mthod(branchtype);
        method=mth_override(method,inp);
    catch
        warning('replace_method_pars:override',...
            'replace_method_pars: overriding method fields for %s failed',branchtype);
        method=inp;
    end
else
    method=inp;
end
if isempty(varargin)
    return
end
fnames=fieldnames(method);
arg=varargin;
for i=1:length(fnames)
    fd=fnames{i};
    [sel,arg]=dde_options_filter(fd,arg);
    method.(fd)=dde_set_options(method.(fd),sel,'pass_on');
end
for i=1:length(fnames)
    fd=fnames{i};
    [method.(fd),arg]=dde_set_options(method.(fd),arg,'pass_on');
end
%% adjust method.stability if discretization changes for stst types
if isfield(method,'stability')
    [sel,other]=dde_options_filter('stability',varargin);
    opt1=dde_set_options(method.stability,sel,'pass_on');    
    options=dde_set_options(opt1,other,'pass_on');
    if isfield(options,'discretization') && ...
            (~isfield(inp.stability,'discretization') ||...
            ~strcmp(options.discretization,inp.stability.discretization))
        tmp=df_mthod('stst',options.discretization);
        [sel,other]=dde_options_filter('stability',varargin);
        method.stability=dde_set_options(tmp.stability,sel,'pass_on');
        method.stability=dde_set_options(method.stability,other,'pass_on');
    end
end
end
%% override defaults
function mth=mth_override(mth,old)
fn=fieldnames(mth);
for i=1:length(fn)
    if isempty(mth) 
        if ~isempty(old) 
            mth=old;
        end
    elseif  length(mth)==length(old) && isstruct(mth(1).(fn{i})) && isfield(old(1),fn{i})
        for k=1:length(mth)
            mth(k).(fn{i})=mth_override(mth(k).(fn{i}),old(k).(fn{i}));
        end
    elseif isfield(old(1),fn{i}) && length(mth)==length(old)
        for k=1:length(mth)
            mth(k).(fn{i})=old(k).(fn{i});
        end
    end
end
end
