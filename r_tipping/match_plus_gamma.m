%% add matching conditions between u_+ and gamma
function prob=match_plus_gamma(prob,gluename,iv,uidx,maps)
prob=coco_add_func(prob,gluename,@glue,@dglue,iv,'zero',...
    'uidx',[uidx.u_plus(maps.u_plus.x1_idx);uidx.u_gamma(maps.u_gamma.x1_idx)]);
end
%%
function [iv,y]=glue(~,iv,u)
dim=length(fieldnames(iv));
[up,ug]=deal(u(1:dim),u(dim+(1:dim)));
y=[up(iv.x)-ug(iv.x);up(iv.p)*ug(iv.q)-up(iv.q)*ug(iv.p)];
end
%%
function [iv,J]=dglue(~,iv,u)
dim=length(fieldnames(iv));
[up,ug]=deal(u(1:dim),u(dim+(1:dim)));
ip=iv;
[ig.x,    ig.p,    ig.q]=deal(...
 ip.x+dim,ip.p+dim,ip.q+dim);
J=zeros(2,2*dim);
J(1,[ip.x,ig.x])=[1,-1];
J(2,[ip.p,    ig.q,     ip.q,     ig.p])=...
 [ug(iv.q),up(iv.p),-ug(iv.p),-up(iv.q)];
end
