function out = coll2msh(msh, tbptol)
if nargin<2
    tbptol=1e-14;
end
if isnumeric(msh)
    tbp=msh;
    ntst=sum(abs(tbp(2:end)-tbp(1:end-1))<tbptol)+1;
    out=reshape(tbp,[],ntst);
else
    out=reshape(msh.mesh.tbp,[],msh.maps.NTST);
end
end