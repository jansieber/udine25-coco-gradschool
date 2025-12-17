function y = slope(data, xbp, T0, T, ~, t)

x_shp= data.coll_seg.maps.xbp_shp; % shape of solutino vector
tbp  = data.coll_seg.mesh.tbp; % base points of mesh
xvals=reshape(xbp,x_shp);
xm=coll_mesh_mat(tbp,t/T)*xvals.'; % obtain interpolation value
xtm=coll_mesh_mat(tbp,t/T,'diff',1)*xvals.';
y=[t;xm(data.idx);xtm(data.idx)];
end
