function [ind,len,xfields]=dde_ind_from_point(point,free_par_ind,xind,base)
%% Determine indices for entries in vector from structure
%
% $Id: dde_ind_from_point.m 336 2019-05-09 17:04:13Z jansieber $
%%
if nargin<2
    free_par_ind=1:size(point.parameter,2);
end
if nargin<4
    base=0;
end
[fields,xfields]=feval(['dde_',point.kind,'_vars']);
if isstruct(fields)
    ind=fields;
else
    [ind,len]=fields{2}(point,free_par_ind,base,xind);
    return
end
if nargin<3
    xind=1:size(point.(xfields{1}),1);
end
point.parameter=1:length(free_par_ind);
for i=1:length(xfields)
    tmp=point.(xfields{i});
    point.(xfields{i})=tmp(xind,:);
end
prev=base;
fnames=fieldnames(fields);
for i=1:length(fnames)
    fd=fnames{i};
    iar=reshape(1:numel(point.(fd)),size(point.(fd)));
    ind.(fd)=prev(end)+iar;
    if ~isempty(iar)
        prev=ind.(fd);
    end
    if fields.(fd)==2
        prev=prev(end)+iar;
        ind.(fd)=struct('re',ind.(fd),'im',prev);
    end
end
len=prev(end)-base;
end