function [x,ind]=dde_x_from_point(points,free_par_ind,xind)
%% extract variable values from point structure
%
% $Id: dde_x_from_point.m 374 2019-09-14 14:02:58Z jansieber $
%%
x=[];
ind=repmat(struct(),1,0);
if isempty(points)
    return
end
pt=points(1);
[ind,len,xfields]=dde_ind_from_point(pt,free_par_ind);
if nargin<2
    free_par_ind=1:size(pt.parameter,2);
end
if nargin<3
    xind=1:size(pt.(xfields{1}),1);
end
x=NaN(len,numel(points));
fnames=fieldnames(ind);
for i=1:length(points)
    pt=points(i);
    pt.parameter=pt.parameter(free_par_ind);
    for j=1:length(xfields)
        tmp=pt.(xfields{j});
        pt.(xfields{j})=tmp(xind,:);
    end
    for k=1:length(fnames)
        fd=ind.(fnames{k});
        if ~isstruct(fd)
            x(reshape(fd,[],1),i)=reshape(pt.(fnames{k}),[],1);
        else
            x(reshape(fd.re,[],1),i)=reshape(real(pt.(fnames{k})),[],1);
            x(reshape(fd.im,[],1),i)=reshape(imag(pt.(fnames{k})),[],1);
        end
    end 
end
x=reshape(x,[len,size(points)]);
end
