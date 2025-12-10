%% add matching conditions between u_+ and gamma
function prob=match_plus_gamma(prob,gluename,uidx,maps,adj)
prob=coco_add_func(prob,gluename,@glue,@dglue,[],'zero',...
    'uidx',[uidx.u_plus(maps.u_plus.x1_idx);uidx.u_gamma(maps.u_gamma.x1_idx)]);
if nargin<5 ||~adj.do
    return
end
aidxrg=[adj.idx.u_plus(adj.opt.u_plus.x1_idx);...
        adj.idx.u_gamma(adj.opt.u_gamma.x1_idx)];
prob=loc_add_adjt(prob,gluename,adj,'aidx',aidxrg);
end
%%
function [d,y]=glue(~,d,u)
[up,ug]=deal(u(1:3),u(4:6));
y=[up(1)-ug(1);up(2)*ug(3)-up(3)*ug(2)];
end
%%
function [d,J]=dglue(~,d,u)
[up,ug]=deal(u(1:3),u(4:6));
J=[1,  0,     0,   -1,   0,    0;...
   0,ug(3),-ug(2),  0,-up(3),up(2)];
end
%%
function [d,J]=d2glue(~,d,u)
J=zeros(2,length(u),length(u));
[J(2,2,6),J(2,6,2)]=deal(1);
[J(2,3,5),J(2,5,3)]=deal(-1);
end
