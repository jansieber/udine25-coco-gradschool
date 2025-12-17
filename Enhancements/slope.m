function y = slope(data, xbp, T0, T, ~, v)

x_shp= data.coll_seg.maps.xbp_shp; % shape of solutino vector
tbp  = data.coll_seg.mesh.tbp; % base points of mesh
xvals=reshape(xbp,x_shp);
xm=coll_mesh_mat(tbp,v/T)*xvals.'; % obtain interpolation value
xtm=coll_mesh_mat(tbp,v/T,'diff',1)*xvals.';
y=[v;xm(data.idx);xtm(data.idx)];
end
