%% general r.h.s. generator
function fcoco=f2coco(f)
fcoco=@(p,d,u)fcoco_gen(p,d,u,f);
end
function [data,y]=fcoco_gen(~,data,u,f)
y=f(u);
end
