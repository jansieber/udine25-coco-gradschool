function tri=triangulate_adjacency(adj,varargin)
default={'level',2};
options=sco_set_options(default,varargin,'pass_on');
np=size(adj,1);
tri=NaN(0,options.level+1);
for i=1:np
    itri=find(adj(i,:));
    if length(itri)<=1
        continue
    end
    nb=adj(itri,itri);
    if options.level==2
        [ir,ic]=find(nb);
        ipairs=unique(sort([itri(ir);itri(ic)],1)','rows');
        tri=cat(1,tri,[repmat(i,size(ipairs,1),1),ipairs]);
    else
        trirec=triangulate_adjacency(nb,options.level-1);
        triloc=cat(2,repmat(i,size(trirec,1),1),itri(trirec));
        tri=cat(1,tri,triloc);
    end
end
tri=unique(sort(tri,2),'rows');
end
