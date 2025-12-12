
function prob = coll_construct_measure(prob, collfid,name)
%% COLL_CONSTRUCT_MEASURES   Add monitoring of solution measures.
% copied and pasted from coll_construct_err.m and coll_add.m->bddat
data=coco_get_func_data(prob,collfid,'data');
seg  = data.coll_seg;
uidx = coco_get_func_data(prob, collfid, 'uidx');
efid = coco_get_id(collfid, {'measures',name});
prob = coco_add_func(prob, efid{1}, @coll_meas, data, ...
  'regular', efid{2}, 'uidx', uidx(seg.maps.xbp_idx), ...
  'remesh', @coll_meas_remesh, 'passChart');
end
%%
function [data, chart, y] = coll_meas(prob, data, chart, u) %#ok<INUSL>
pr = data.pr;
seg = pr.coll_seg;
fid = seg.fid;
cdata = coco_get_chart_data(chart, fid);
if isfield(cdata, 'measures')
  y = cdata.measures;
else
  maps = seg.maps;
  mesh  = seg.mesh;
  x     = u;
  xcn   = reshape(maps.W*x, maps.x_shp);
  NTST  = maps.NTST;
  wts   = mesh.gwt; % Gauss weights
  kas   = mesh.gka; % Warping coefficients
  nrmx  = sqrt((0.5/NTST)*sum(xcn.*xcn,1)*(wts.*kas)');
  y  = nrmx;
  cdata.measures = y;
  chart = coco_set_chart_data(chart, fid, cdata);
end
end
%%
function [prob, stat, xtr] = coll_meas_remesh(prob, data, chart, ub, Vb) %#ok<INUSD>
seg  = data.coll_seg;
maps = seg.maps;
xtr  = [];
uidx = coco_get_func_data(prob, seg.fid, 'uidx');
prob = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
stat = 'success';
end
