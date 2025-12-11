%% add matching conditions between u_+ and gamma
function prob=match_plus_gamma(prob,gluename,uidx,maps)
prob=coco_add_func(prob,gluename,@glue,@dglue,[],'zero',...
    'uidx',[uidx.u_plus(maps.u_plus.x1_idx);uidx.u_gamma(maps.u_gamma.x1_idx)]);
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
