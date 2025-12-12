function [adj,id,onesided]=adjacency_from_atlas(input,varargin)
%% create (full) adjacency matrix from neighborhood information of atlas
% input can be run name, atlas, or neighborhood cell array
default={'onesided','add'};
options=sco_set_options(default,varargin,'pass_on');
atlas=[];
if isstruct(input)
    atlas=input;
end
if ischar(input)
    atlas=coco_bd_read(input,'atlas');
end
if ~isempty(atlas)
    nbh=cellfun(@(x)x.neigh,atlas.charts,'uniformoutput',false);
    id=cellfun(@(x)x.id,atlas.charts);
    nbh=nbh(id);
else
    nbh=input;
end
adj=zeros(length(nbh));
for i=1:length(nbh)
    for k=1:length(nbh{i})
        if nbh{i}(k)==0
            continue
        end
        adj(i,nbh{i}(k))=adj(i,nbh{i}(k))+1;
    end
end
[irs,ics]=find(adj~=adj');
onesided=reshape([irs,ics],[],2);
switch options.onesided
    case 'remove'
        adj(irs,ics)=0;
    case 'add'
        adj_add=sparse(irs,ics,ones(size(irs,1),1),size(adj,1),size(adj,2));
        adj=max(adj,adj_add);
end
onesided=unique(sort(onesided,2),'rows');
end
