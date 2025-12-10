%% add condition relating T to r
function [prob,rtini]=fix_r_T(prob,rtname,uidx,u0,maps,ip,adj)
segs={'u_plus','u_minus'};
d0=struct('rt0',0);
for i=1:length(segs)
    u0i=u0.(segs{i});    
    map=maps.(segs{i});
    ind=uidx.(segs{i});
    idx=[map.p_idx(ip.r);map.T_idx];
    [~,rt0]=r_T(prob,d0,u0i(idx));
    rtini.(segs{i})=rt0;
    d=struct('rt0',rt0);
    id=coco_get_id(rtname,segs{i});
    prob=coco_add_func(prob,id,@r_T,d,'zero','uidx',ind(idx),'f+df');
    if nargin>6 && adj.do
        aind=adj.idx.(segs{i});
        amap=adj.opt.(segs{i});
        idx=[amap.p_idx(ip.r);amap.T_idx];
        prob=loc_add_adjt(prob,id,adj,'aidx',aind(idx),'f+df');
    end
end
end
%%
function [d,y,J]=r_T(~,d,u)
y=u(1)*u(2)-d.rt0;
if nargout<3
    return
end
J=[u(2),u(1)];
end

