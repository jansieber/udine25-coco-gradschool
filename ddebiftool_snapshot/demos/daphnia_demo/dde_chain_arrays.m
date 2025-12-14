function [pt,istruc]=dde_chain_arrays(tab,name,point)
fun=tab.fun_ref(tab.fun_names.(name));
nt=size(point.profile,2);
point.parameter=point.parameter(fun.par);
point.profile=point.profile(fun.x,:);
[tau,xxd]=dde_coll_delays(fun.funcs,point,'t',point.mesh,'max_xderiv',dde_num_delays(fun.funcs,'max_xderiv'));
parg=repmat(point.parameter,1,1,nt);
[r,yarr,istruc]=fun.funcs.wrap_rhs(xxd,parg);
hst=istruc.histories;
hnames=fieldnames(hst.names);
[taumesh,timemesh]=ndgrid(istruc.msh.t,point.mesh*point.period);
ibd=istruc.args.ibd;
switch ibd.type{1}{1}
    case 'parameter'
        taumesh=taumesh*point.parameter(ibd.local.parameter{1});
    case 'x'
        taumesh=taumesh.*repmat(point.profile(ibd.local.x{1},:),size(taumesh,1),1);
end
pt=dde_coll_create('mesh',struct('time',timemesh,'tau',taumesh),'degree',istruc.msh.degree);
ys=struct();
for i=1:length(hnames)
    ys.(hnames{i}).int=yarr(hst.xvals{1+hst.names.(hnames{i})}{hst.int},:,:);
    ys.(hnames{i}).id=yarr(hst.xvals{1+hst.names.(hnames{i})}{hst.id},:,:);
    if length(hst.xvals{1+hst.names.(hnames{i})}{hst.int})==1
        ys.(hnames{i}).int=reshape(ys.(hnames{i}).int,[],nt);
    end
    if length(hst.xvals{1+hst.names.(hnames{i})}{hst.id})==1
        ys.(hnames{i}).id=reshape(ys.(hnames{i}).id,[],nt);
    end
end
pt.profile=ys;
end