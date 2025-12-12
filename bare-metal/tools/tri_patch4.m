function [fixed,adj,info]=tri_patch4(adj,varargin)
tri=triangulate_adjacency(adj);
tri=sort(tri,2);
trc=tri_count(tri);
%% check for redundant triangles
nadj=max(tri(:));
sp=@(irc)sparse(irc(:,1),irc(:,2),ones(size(irc,1),1),nadj,nadj);
spsym=@(a)a+a';
sps=@(irc)spsym(sp(irc));
adj3=sp(trc{3});
info.redundant_tri=triangulate_adjacency(adj3,'level',2);
info.redundant_edges=list_edges(info.redundant_tri);
info.flaws=setdiff(trc{3},info.redundant_edges,'rows');
info.flaws=cat(1,info.flaws,trc{4:end});
frnodes=unique(trc{1}(:));
bd_ind=1:length(frnodes);
mfr=zeros(1,max(frnodes));
mfr(frnodes)=bd_ind;
ifr=mfr(trc{1});
ifr2=cat(1,ifr,ifr(:,[2,1]));
adj_fr=sparse(ifr2(:,1),ifr2(:,2),ones(size(ifr2,1),1));
g=graph(adj_fr);
cc4=find_comps(g,4);
crossings=find(g.degree()>2);
assert(all(mod(g.degree(crossings),2)==0));
gcut=rmnode(g,crossings);
cc3=find_comps(gcut,3);
cc=[cc4,cc3];
indc=cellfun(@find_dist,cc,'uniformoutput',false);
ind=reshape(cat(1,indc{:}),[],2);
info.addlinks=reshape([frnodes(ind(:,1)),frnodes(ind(:,2))],[],2);
adj_sub=sps(info.redundant_edges);
adj_add=sps(info.addlinks);
adj=max(max(adj-adj_sub,0),adj_add);
fixed=triangulate_adjacency(adj);
fixed=sortrows(sort(fixed,2));
end
function gc=find_comps(gr,sz)
ccnum=gr.conncomp();
cc=cell(max(ccnum),1);
for i=1:max(ccnum)
    cc{i}=find(ccnum==i);
end
loops=find(cellfun(@(c)length(c)==sz,cc));
cc=cc(loops);
gc=cell(1,0);
for i=length(cc):-1:1
    gc{i}={subgraph(gr,cc{i}),cc{i}};
end
end
function irc=find_dist(grc)
[ir,ic]=find(grc{1}.distances==2,1,'first');
if ~isempty(ir)
    irc=grc{2}([ir,ic]);
else
    irc=zeros(0,2);
end
end
function trc=tri_count(tri)
%% count in how many triangles each edge is involved
trimat=list_edges(tri);
[~,ia,ix]=unique(trimat,'rows');
iad=diff([ia;size(trimat,1)+1].');
for i=max(iad):-1:1
    trc{i}=trimat(ia(iad==i),:);
end
%tr=triangulation(tri,xp(:,1),xp(:,2),xp(:,3));
%fr=tr.freeBoundary();
end
function trimat=list_edges(tri)
trimat=permute(reshape(tri(:,[1,2,2,3,3,1]),size(tri,1),2,[]),[1,3,2]);
trimat=sortrows(sort(reshape(trimat,[],2),2));
end