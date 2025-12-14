function br_warn(method,conds,str,varargin)
%% warning messages during continuation depending on flags
if ~isfield(method,'warnings') || ~method.warnings
    return
end
for i=1:length(conds)
    if ischar(conds{i}) && ~isfield(method,conds{i})
        return
    elseif islogical(conds{i}) && ~conds{i}
        return
    end
end
fprintf(str,varargin{:});
end