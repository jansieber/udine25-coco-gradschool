function point=dde_TorusBifurcation_postprocess(point,data)
%% renormalize eigenvector to be as close as possible to reference
%%
prev=data.previous;
n=size(prev.profile,1)/3;
kn=@(mat)kron(mat,eye(n));
wabs=kn(diag([0,1,1]));
wph=kn(blkdiag(0,[0,1;-1,0]));
sc=@(p1,p2,mat)dde_coll_profile_dot(p1,p2,'inner_matrix',mat);
pnrm=sc(point,point,wabs);
point.profile=kn(diag([1,1/pnrm,1/pnrm]))*point.profile;
c=sc(point,prev,wabs)+1i*sc(point,prev,wph);
phi=angle(c);
krot=kn(blkdiag(1,[cos(phi),-sin(phi);sin(phi),cos(phi)]));
point.profile=krot*point.profile;
end

