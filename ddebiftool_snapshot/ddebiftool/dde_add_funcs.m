function tab=dde_add_funcs(tab,name,funcs,varargin)
%% add new equations and delays to tab (labelled by name)
% args x_old, par_old, rhstau_old are index sets pointing to previously
% added tab entries rhstau is number of new columns of xx in sys_rhs and
% sys_tau, delays are also numbered in the order they come in. nf is number
% of equations, x and par are values for new variables that sys_rhs and
% sys_tau depend on, determining dimensions and new indices by their size.
default={...
    'x',zeros(0,1),'par',zeros(1,0),...
    'x_old',zeros(0,1),'par_old',zeros(1,0),...
    'rhstau_args_old',zeros(0,1),'rhstau_args',dde_num_delays(funcs),...
    'nf',NaN,'checkrhs',true};
options=dde_set_options(default,varargin,'pass_on');
tab=dde_funtab(tab);
xtest=cat(1,tab.xtest(options.x_old),options.x);
partest=cat(2,tab.partest(options.par_old),options.par);
x_ind=[options.x_old(:)',length(tab.xtest)+(1:length(options.x))];
par_ind=[options.par_old(:)',length(tab.partest)+(1:length(options.par))];
rhstau_args_ind=[options.rhstau_args_old(:)',tab.ntau+(1:options.rhstau_args)];
tauprov=tab.ntau+(1:dde_num_delays(funcs));
nf=get_nf(funcs,xtest,partest,options);
f_ind=tab.nf+(1:nf);
tab.xtest=cat(1,tab.xtest,options.x);
tab.partest=cat(2,tab.partest,options.par);
funcs.xpattern=dde_funcs_xpattern(funcs,length(x_ind));
fun=dde_funformat_create('funcs',funcs,'par',par_ind,'x',x_ind,'rhstau_args',rhstau_args_ind,...
    'f',f_ind,'tauprovided',tauprov);
fun_num=length(tab.fun_ref)+1;
tab.fun_ref(fun_num)=fun;
tab.fun_names.(name)=fun_num;
tab.fun_namelist(end+1)={name};
tab.ntau=tab.ntau+length(tauprov);
tab.taumap(tauprov,1)=fun_num;
tab.nf=tab.nf+nf;
end

function nf= get_nf(funcs,xtest,partest,options)
assert(options.checkrhs || ~isnan(options.nf));
nf=options.nf;
if options.checkrhs
    [ntaum1,max_xderiv]=dde_num_delays(funcs);
    xtest=xtest(:,ones(ntaum1+1,1),ones(max_xderiv+1,1));
    y=funcs.wrap_rhs(xtest,partest);
    nf=size(y,1);
end
if isnan(options.nf)
    options.nf=nf;
end
assert(options.nf==nf);
end