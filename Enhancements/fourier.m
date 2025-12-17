function y = fourier(data, xbp, ~, ~, ~)
%% compute Fourier modes of xbp
mps = data.coll_seg.maps;
msh = data.coll_seg.mesh;
% matrix corresponding to linear map x->int_0^1 x(s)ds
intJ=[-1,1]*coll_mesh_int(msh.tbp,[0,1]);
% turn into diagonal matric of intergration weights
Wt=diag(sparse(intJ));
% reshape xbp to form nx x ntst*(deg+1)
xbp=reshape(xbp,mps.xbp_shp);
ebp=exp(-2i*pi*msh.tbp*data.n);
y=xbp*Wt*ebp;
y=y(:);
end
