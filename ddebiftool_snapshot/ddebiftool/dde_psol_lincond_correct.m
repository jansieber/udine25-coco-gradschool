%% correct projections and transformations for psol bifurcations
function S=dde_psol_lincond_correct(p,Sin)
S=Sin;
profiledim=size(p.profile,1);
[nc,xdim]=size(S.stateproj);
nv=profiledim/xdim-1;
switch S.fieldname
    case 'x'
        S.stateproj=cat(2,S.stateproj,zeros(nc,profiledim-xdim));
        S.fieldname='profile';
    case 'v'
        S.stateproj=kron(eye(nv),S.stateproj);
        S.stateproj=cat(2,zeros(nc*nv,xdim),S.stateproj);
        S.trafo=kron(eye(nv),S.trafo);
        S.condprojmat=kron(eye(nv),S.condprojmat);
        S.fieldname='profile';
end
end
