function [tri,adj]=triangles_from_atlas(varargin)
%% create trimesh from atlas
default={'remove_redundant',true,'run',[],'atlas',[]};
options=sco_set_options(default,varargin,'pass_on');
% input can be run name or atlas
if isempty(options.atlas)
    assert(~isempty(options.run),'triangles_from_atlas:arguments',...
        'triangles_from_atlas: provide one of arguments run or atlas');
    options.atlas=coco_bd_read(options.run,'atlas');
end
nbh=cellfun(@(x)x.neigh,options.atlas.charts,'uniformoutput',false);
adj=adjacency_from_atlas(nbh);
%% if requested, remove redundant triangles
% these are those that involve points which have more than 2 neighbors and
% for which all neighbors are mutually connected
if options.remove_redundant
for i=1:length(nbh)
    itri=find(adj(i,:));
    if isempty(itri)
        continue
    end
    nb=adj(itri,itri);
    if length(itri)>2 && min(min(nb+diag(ones(size(nb,1),1))))>0
        adj(i,:)=0;
        adj(:,i)=0;
    end
end
%% extract triangles
tri=NaN(0,3);
for i=1:length(nbh)
    itri=find(adj(i,:));
    if isempty(itri)
        continue
    end
    nb=adj(itri,itri);
    [ir,ic]=find(nb);
    ipairs=unique(sort([itri(ir);itri(ic)],1)','rows');
    tri=cat(1,tri,[repmat(i,size(ipairs,1),1),ipairs]);
end
tri=unique(sort(tri,2),'rows');
end

