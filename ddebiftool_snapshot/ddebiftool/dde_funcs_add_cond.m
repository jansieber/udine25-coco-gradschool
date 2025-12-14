function funcs=dde_funcs_add_cond(funcs,sys_cond,varargin)
default={'reference',false,'name','sys_cond'};
options=dde_set_options(default,varargin,'pass_on');
if ~isstruct(sys_cond)
    if ~iscell(sys_cond)
        sys_cond={sys_cond};
    end
    sys_cond_struc=repmat(dde_sys_cond_create(),0,1);
    for i=length(sys_cond):-1:1
        sys_cond_struc(i,1)=dde_sys_cond_create('name',options.name,...
            'fun',sys_cond{i},'reference',options.reference);
    end
else
    sys_cond_struc=sys_cond(:);
end
if ~isstruct(funcs.sys_cond)
    funcs.sys_cond=dde_sys_cond_create('name',options.name,'fun',funcs.sys_cond,...
        'reference',funcs.sys_cond_reference);
end
funcs.sys_cond=cat(1,funcs.sys_cond,sys_cond_struc);
end