function plot_chart(c,ind,clrnum)
x0=c.xp(ind);
vd=reshape([c.P.v{:}],c.dim,[]);
nv=size(vd,2);
pvert=c.xp+c.TSp*vd;
for k=c.P.nFaces:-1:1
    pfacemid(:,k)=mean(pvert(:,c.P.faceV{k}),2);
end
pfacenbh=facenbh(c.P.faceV,c.dim);
pvert=pvert(ind,:);
nv=size(pfacemid,2);
pfacemid=pfacemid(ind,:);
chain=find_chain(pfacenbh);
chain=[chain,chain(1)];
ax=gca;
clr=lines();
if nargin<3
    clrnum=get(ax,'LineStyleOrderIndex');
end
vs=sortvertices(c.P);
vs=[vs(:)',vs(1)];
patch(ax,'Faces',1:nv,'Vertices',pvert(:,vs)',...
    'FaceColor', clr(clrnum,:), 'FaceAlpha', 0.5);
ish=ishold(ax);
hold(ax,'on');
plot3(ax,pfacemid(1,chain),pfacemid(2,chain),pfacemid(3,chain),'x-','linewidth',2);
plot3(ax,pvert(1,vs),pvert(2,vs),pvert(3,vs),'o-','linewidth',2);
clrnum=get(ax,'LineStyleOrderIndex');
plot3(ax,x0(1),x0(2),x0(3),'o','color','k','linewidth',2,...
    'MarkerFaceColor',clr(clrnum,:),'MarkerSize',8);
if ~ish
    hold(ax,'off');
end
end
function vs = sortvertices(P)

temp = cell2mat(P.faceV');
vs   = temp(1,:);
lookfor   = temp(1,2);
temp(1,:) = [0, 0];
for k=2:P.nFaces-1
  i = find(lookfor==temp(:,1));
  if isempty(i)
    i = find(lookfor==temp(:,2));
    lookfor = temp(i(1),1);
  else
    lookfor = temp(i(1),2);
  end
  vs(k+1) = lookfor;
  temp(i(1),:) = [0, 0];
end

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