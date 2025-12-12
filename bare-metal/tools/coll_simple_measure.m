function prob = coll_simple_measure(prob, collfid,name)
%% COLL_CONSTRUCT_MEASURES   Add monitoring of solution measures.
% copied and pasted from coll_construct_err.m and coll_add.m->bddat
efid = coco_get_id(collfid, {'measures',name});
colldata=coco_get_func_data(prob,collfid,'data');
seg  = colldata.coll_seg;
uidx = coco_get_func_data(prob, collfid, 'uidx');
data=struct('fid',collfid);
prob = coco_add_func(prob, efid{1}, @coll_meas, data, ...
  'regular', efid{2}, 'uidx', uidx(seg.maps.xbp_idx), ...
  'remesh', @coll_meas_remesh);
end
%%
function [data, y] = coll_meas(prob, data, u)
colldata=coco_get_func_data(prob,data.fid,'data');
seg = colldata.pr.coll_seg;
maps = seg.maps;
mesh  = seg.mesh;
x     = u;
xcn   = reshape(maps.W*x, maps.x_shp);
NTST  = maps.NTST;
wts   = mesh.gwt; % Gauss weights
kas   = mesh.gka; % Warping coefficients
nrmx  = sqrt((0.5/NTST)*sum(xcn.*xcn,1)*(wts.*kas)');
y  = nrmx;
end
function [prob, stat, xtr] = coll_meas_remesh(prob, data, chart, ub, Vb) %#ok<INUSD>
colldata=coco_get_func_data(prob,data.fid,'data');
seg  = colldata.coll_seg;
maps = seg.maps;
xtr  = [];
uidx = coco_get_func_data(prob, seg.fid, 'uidx');
prob = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
stat = 'success';
end
