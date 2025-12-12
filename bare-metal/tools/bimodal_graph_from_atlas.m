function gr=bimodal_graph_from_atlas(atlas,info,varargin)
default={'maxdefect',Inf};
options=sco_set_options(default,varargin,'pass_on');
%% find, truncate and symmetrize adjacency matrix
adj=adjacency_from_atlas(info.nbh);
adj=double(adj>0&adj'>0);
adj=triu(adj,1);
np=info.npts;
%% for each point determine polyhedron vertices (from plot_atlas_kd)
% and midpoints of faces, also determine, which faces are neighbors
pvert=cell(1,np);
pfacenbh=cell(1,np);
pfacemid=cell(1,np);
pface2index=cell(1,np);
index2pface=NaN(sum(cellfun(@(c)c.P.nFaces,atlas.charts))+np,2);
index2pface(1:np,1)=1:np;
index=np;
for i=1:np
    c=atlas.charts{i};
    vd=reshape([c.P.v{:}],info.dim,[]);
    nv=size(vd,2);
    pvert{i}=info.xp(i*ones(1,nv),:)'+c.TSp*vd;
    for k=c.P.nFaces:-1:1
        pfacemid{i}(:,k)=mean(pvert{i}(:,c.P.faceV{k}),2);
    end
    pfacenbh{i}=facenbh(c.P.faceV,info.dim);
    index2pface(index+(1:c.P.nFaces),1)=i;
    index2pface(index+(1:c.P.nFaces),2)=1:c.P.nFaces;
    pface2index{i}=index+(1:c.P.nFaces);
    index=index+c.P.nFaces;
end
index2pface=index2pface(1:index,:);
%% for each neighboring pair associate corresponding face mid points
% for identification
[ir,ic]=find(adj);
%identify=NaN(2,2,length(ir));
defect=NaN(1,length(ir));
newfacemid=pfacemid;
newindex2face=[index2pface,NaN(size(index2pface))];
%rejected=false(length(ir),1);
for i=1:length(ir)
    n1=atlas.charts{ir(i)}.neigh;
    n2=atlas.charts{ic(i)}.neigh;
    i1=find(n1==ic(i));
    i2=find(n2==ir(i));
    %pfm1=pfacemid{ir(i)};
    %pfm2=pfacemid{ic(i)};
    %cost=sco_match_vector(pfm1,pfm2,'costmatonly',true);
    %[~,ilin]=min(cost(:));
    %[i1,i2]=ind2sub(size(cost),ilin);
    %identify(:,:,i)=[ir(i),i1;ic(i),i2];
    defect(i)=norm(pfacemid{ic(i)}(:,i2)-pfacemid{ir(i)}(:,i1));
    if defect(i)<=options.maxdefect
        newmid=mean([pfacemid{ic(i)}(:,i2),pfacemid{ir(i)}(:,i1)],2);
        newfacemid{ir(i)}(:,i1)=newmid;
        newfacemid{ic(i)}(:,i2)=newmid;
        newindex2face(pface2index{ic(i)}(i2),:)=NaN;
        newindex2face(pface2index{ir(i)}(i1),3:4)=[ic(i),i2];
    %else
    %    rejected(i)=true;
    end
end
%identify=identify(:,:,defect<=options.maxdefect);
newindex2face=newindex2face(~all(isnan(newindex2face),2),:);
newface2index=cell(1,np);
index2face=cell(size(newindex2face,1),1);
for i=1:size(newindex2face,1)
    i2f=cell(1,0);
    if ~isnan(newindex2face(i,2))
        newface2index{newindex2face(i,1)}(newindex2face(i,2))=i;
        i2f{1}=[newindex2face(i,1);newindex2face(i,2)];
    end
    if ~isnan(newindex2face(i,3))
        newface2index{newindex2face(i,3)}(newindex2face(i,4))=i;
        i2f{2}=[newindex2face(i,3);newindex2face(i,4)];
    end
    index2face{i}=i2f;
end
%% assemble graph connecting basepoints to faces and faces to each other
n_be=length(index2face);
nlinks=sum(cellfun(@(p)size(p,1),pfacenbh))+sum(cellfun(@length,newfacemid));
irc=NaN(nlinks,2);
index=0;
for i=1:np
    nnb=size(pfacenbh{i},1);
    irc(index+(1:nnb),:)=newface2index{i}(pfacenbh{i});
    index=index+nnb;
    nf=size(newfacemid{i},2);
    irc(index+(1:nf),:)=[i*ones(nf,1),newface2index{i}'];
    index=index+nf;
end
irc=irc(1:index,:);
x=NaN(size(info.xp,2),length(index2face));
for i=1:length(index2face)
    i2f=index2face{i};
    if isempty(i2f)
        x(:,i)=info.xp(i,:)';
    else
        i2f=i2f{1};
        x(:,i)=newfacemid{i2f(1)}(:,i2f(2));
    end
end
%% create adjacency matrix for bimodal graph
b_adj=sparse(irc(:,1),irc(:,2),ones(size(irc,1),1),n_be,n_be);
b_adj=b_adj+b_adj';
gr=struct(...
    'adjacency',b_adj,...
    'face2index',{newface2index},...
    'index2face',{index2face},...
    'face_neighbors',{pfacenbh},...
    'xpoints',x);
end
function nbh=facenbh(faceV,dim)
if isempty(faceV)
    nbh=NaN(0,2);
    return
end
nbhar=false(length(faceV),length(faceV));
for i=1:length(faceV)
    for j=i+1:length(faceV)
        if length(intersect(faceV{i},faceV{j}))>=dim-1
            nbhar(i,j)=true;
        end
    end
end
[nbh(:,1),nbh(:,2)]=find(nbhar);
end
