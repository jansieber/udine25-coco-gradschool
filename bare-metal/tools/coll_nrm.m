function y = coll_nrm(data, xbp, T0, T, p, t)
seg = data.coll_seg;
maps = seg.maps;
mesh  = seg.mesh;
x     = xbp;
xcn   = reshape(maps.W*x, maps.x_shp);
NTST  = maps.NTST;
wts   = mesh.gwt; % Gauss weights
kas   = mesh.gka; % Warping coefficients
nrmx  = sqrt((0.5/NTST)*sum(xcn.*xcn,1)*(wts.*kas)');
y  = nrmx;
end