function pofoldini=dde_pofold_template(funcs,point)
pofoldini=point;
ip=funcs.ip;
pofoldini.profile(ip.dim+(1:ip.dim),:)=0;
pofoldini.parameter([ip.beta,ip.nullparind(:,2)'])=0;
end