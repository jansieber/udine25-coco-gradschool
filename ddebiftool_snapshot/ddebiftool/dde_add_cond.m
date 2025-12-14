function funcs=dde_add_cond(name,funcs,options,testpoint,condname,output)
if nargin<5
    condname='usercond';
end
if nargin<6
    output=true;
end
%% add new sys_cond if given by user
if ~isempty(options.(condname))
    if output && ~options.outputfuncs
        error([name,':sys_cond'],[name,': ',...
            'new extra conditions change functions structure. Set\n',...
            '''outputfuncs'' to  true.']);
    end
    options.usercond=dde_test_cond(options.usercond,testpoint,'name',condname);
    funcs=dde_funcs_add_cond(funcs,options.usercond);
end
end
