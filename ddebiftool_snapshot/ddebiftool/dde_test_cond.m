function out=dde_test_cond(cond,testpoint,varargin)
default={'name',''};
options=dde_set_options(default,varargin,'pass_on');
if isstruct(cond) || isempty(cond)
    out=cond;
    return
end
if ~iscell(cond)
    cond={cond};
end
for i=numel(cond):-1:1
    try
        [rdum,Jdum]=cond{i}(testpoint); %#ok<ASGLU>
        reference=false;
    catch
        reference=true;
    end
    if isempty(options.name)
        name=func2str(cond{i});
    else
        name=options.name;
    end
    out(i)=dde_sys_cond_create('name',name,...
        'fun',cond{i},'reference',reference);
end
out=out(:);
end
