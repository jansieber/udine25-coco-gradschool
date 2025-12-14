function [args,dargs]=dde_cell_from_arrays(xx,xdev, par,pdev,...
    nx,ndelays,npar,mfrep,iscollected)
%% convert arrays into cells for dde_sym_rhs_wrap and dde_sym_tau_wrap
if ~iscollected
    args={xx.',par.'};
    dargs={xdev.',pdev.'};
    return
end
for i=npar:-1:1
    parc{i}=par(i,:);
end
for i=npar*mfrep:-1:1
    dparc{i}=pdev(i,:);
end
for i=nx*ndelays:-1:1
    xxc{i}=xx(i,:);
end
for i=nx*ndelays*mfrep:-1:1
    dxxc{i}=xdev(i,:);
end
args=[xxc,parc];
dargs=[dxxc,dparc];
end