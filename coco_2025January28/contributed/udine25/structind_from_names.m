function [s,nnames]=structind_from_names(names)
%% generate structure with counts from cell array of names
nnames=length(names);                            % number of names
indc=[names;num2cell(1:nnames)];
s=struct(indc{:});        % ip.r will be number of parameter r in array pars
end