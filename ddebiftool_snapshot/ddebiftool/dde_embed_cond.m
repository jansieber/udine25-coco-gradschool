function embcond=dde_embed_cond(funcs,origcond,part,varargin)
if ~isfield(funcs,'get_comp')||~isfield(funcs,'embed')
    embcond=repmat(dde_sys_cond_create(),0,1);
    if ~isempty(origcond)
        warning('dde_embed_cond:conversion',['dde_embed_cond: ',...
            'non-empty condition requested for embedding, but get_comp ',...
            'or embed misssing']);
    end
    return
end
default={'reference',false};
options=dde_set_options(default,varargin,'pass_on');
if isfield(funcs,'userfuncs')
    reference=funcs.userfuncs.sys_cond_reference;
else
    reference=options.reference;
end
if ~isstruct(origcond)
    origcond=dde_sys_cond_create('name','sys_cond','fun',origcond,'reference',reference);
end
embcond=origcond;
for i=length(origcond):-1:1
    embfun=@(p,pref)loc_embcond(p,pref,origcond(i),funcs.get_comp,funcs.embed,part);
    embcond(i)=dde_sys_cond_create('name',[part,'.',origcond.name],'fun',embfun,'reference',true);
end
end
%%
function [r,J]=loc_embcond(p,pref,origcond,get_comp,embed,part)
userpoint=get_comp(p,part);
userref=get_comp(pref,part);
if origcond.reference
    args={userpoint,userref};
else
    args={userpoint};
end
[r,Ju]=origcond.fun(args{:});
J=embed(Ju,part,pref);
end