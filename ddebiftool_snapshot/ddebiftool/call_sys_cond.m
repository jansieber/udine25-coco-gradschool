function [r,J,name]=call_sys_cond(funcs,p,pref)
%% call all extra conditions ttached to funcs structure
if nargin<3
    pref=p;
end
if ~isstruct(funcs.sys_cond)
    if isfield(funcs,'sys_cond_reference') && funcs.sys_cond_reference
        args={p,pref};
    else
        args={p};
    end
    [r,J]=funcs.sys_cond(args{:});
    r=r(:);
    J=J(:);
    name={'sys_cond',1:length(r)};
else
    reference{1}={p};
    reference{2}={p,pref};
    count=0;
    ncond=length(funcs.sys_cond);
    [rc,Jc]=deal(cell(ncond,1));
    rc{1}=zeros(0,1);
    Jc{1}=repmat(p,0,1);
    name=cell(ncond,2);
    for i=1:length(funcs.sys_cond)
        args=reference{1+double(funcs.sys_cond(i).reference)};
        [rc{i},Jc{i}]=funcs.sys_cond(i).fun(args{:});
        rc{i}=rc{i}(:);
        Jc{i}=Jc{i}(:);
        name(i,:)={funcs.sys_cond(i).name,count+(1:length(rc{i}))};
        count=count+length(rc{i});
    end
    r=cat(1,rc{:});
    J=cat(1,Jc{:});
end
end