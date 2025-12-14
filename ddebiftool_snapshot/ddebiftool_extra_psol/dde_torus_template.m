function trini=dde_torus_template(funcs,point)
ip=funcs.ip;
trini=point;
trini.profile(ip.null,:)=0;
trini.profile(ip.xrg,:)=point.profile;
trini.parameter(1:ip.nuserpar)=trini.parameter;
trini.parameter(ip.omega)=0.25;
trini.parameter(ip.period)=trini.period;
end