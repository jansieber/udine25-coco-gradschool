function prob = query_construct_save(prob,name)
%% COLL_CONSTRUCT_MEASURES   Add monitoring of solution measures.
% copied and pasted from coll_construct_err.m and coll_add.m->bddat
efid = coco_get_id(name);
prob = coco_add_chart_data(prob, efid, struct(), struct());
data=struct('fid',efid);
prob = coco_add_func(prob, efid, @save_query,data, ...
  'regular',{}, 'uidx',[], 'passChart');
%prob = coco_add_slot(prob, efid{1}, @coco_save_data, data, 'save_full');
end
%%
function [data, chart, y] = save_query(prob, data, chart, u)
cdata = coco_get_chart_data(chart, data.fid);
y = u;
if ~isfield(cdata, 'tab') && isfield(prob.efunc,'acp_idx')
    cdata.idx=struct('idx2par',{prob.efunc.idx2par},...
        'pidx2midx',prob.efunc.pidx2midx,...
        'pidx2fidx',prob.efunc.pidx2fidx,...
        'acp_idx',prob.efunc.acp_idx);
    chart = coco_set_chart_data(chart, data.fid, cdata);
end
end
