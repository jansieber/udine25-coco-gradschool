function tab=dde_close_chain_delay(tab,ind,varargin)
%% add nested distributed delay by adding equation of the form
% h_{j+1}(t-s) = \int_0^s g_j(s,x_j(t-s),par) ds for j=0..k
% 0 = \int_0^(tau*t(i)) h_{j_i}(s) ds - d_i 
%
% where h0=x0=x and xj contains arbitrary combination of h_i with i<j
%
% 0 = int_0^inf g(s,x(t-s),par) exp(-s/tau) (s/tau)^alpha ds - d 
% or equivalently
% 0 = int_0^inf tau g(tau s,x(t-tau s),par) exp(-s) (s)^alpha ds - d 
% (generalized) Laguerre approximation
%
% idist are index of d in variable array, ibd is
% index of tau in variable array, ginp is kernel function g(s,x,p), ntau is
% the delay number after which the new delays (nodes in the integral)
% should be appended (ntau+(1:length(grid.t))
%
% optional argument 'par' passes index of parameters in argument p for g
% (default same as for DDE rhs), 'gt' or 'g_of_t_x flags that first argument t
% of g is present (default true)
default={'name','chain_delay',{'ival','value'},NaN,{'it','itval','t_indices'},NaN,...
    {'t_type'},{'NaN'},...
    {'M','mat','matrices','S'},NaN,{'ix','x_args','x_arguments'},NaN,...
    'maxorder',2};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
ind.args.ival.global={options.ival(:)'};
ind.args.ival.type={repmat({'x'},1,length(options.ival))};
ind.args.t.global={options.it(:).'};
ind.args.t.type={options.t_type(:).'};
ind.args=arg_wrap2cell(ind.args);
[ind.args,ipmap]=local_from_global(ind.args,'parameter');
[ind.args,ixmap]=local_from_global(ind.args,'x');
ind.args.init=init_xpar_reduce(ind.args.init);
ind.sum.input=int_delay_chain_inputs(ind.histories,options.ix);
ind.sum.coeffs_t=prod_deriv_coeffs(options.maxorder*[1,1],options.maxorder);
ind.sum.coeffs_whole=prod_deriv_coeffs([options.maxorder,1],options.maxorder);
%% set matrix by which all delays get combined
ind.sum.M=loc_set_matrix(options.M, ind);
%%
ix0_g=ind_local_x(ind.args.igx);
[ix_g,it_g]=ndgrid(ix0_g(:)',ind.msh.tau);
[ix_bd,ix_init,ix_val,ix_t]=deal(ind_local_x(ind.args.ibd),...
    ind_local_x(ind.args.init),ind_local_x(ind.args.ival),ind_local_x(ind.args.t));
on=@(ix)ones(length(ix),1);
xpattern=dde_join_xpattern([...
    ix_g(:),it_g(:);...
    ix_bd,  on(ix_bd);...
    ix_init,on(ix_init);...
    ix_val, on(ix_val);...
    ix_t,   on(ix_t)]');
%%
tsel=   @(x,p)int_delay_xparchoice(x,p,ind.args.t);
bdsel=  @(x,p)int_delay_xparchoice(x,p,ind.args.ibd);
initsel=@(x,p)int_delay_xparchoice(x,p,ind.args.init);
valsel= @(x)  int_delay_xparchoice(x,zeros(1,0,size(x,3)),ind.args.ival);
nt=ind.args.t.nval;
rhs=@(x,p)int_delay_chain(ind,0,...
    x,p,bdsel(x,p),initsel(x,p),tsel(x,p),valsel(x),...
    zeros(size(x)),zeros(size(p)),zeros(1,size(x,3)),...
    zeros(ind.args.init.nval,size(x,3)),zeros(nt,size(x,3)),...
    zeros(ind.args.ival.nval,size(x,3)));
drhs=@(ord,x,p,dx,dp)int_delay_chain(ind,ord,...
     x, p,bdsel( x, p),initsel( x, p),tsel( x, p),valsel(x),...
    dx,dp,bdsel(dx,dp),initsel(dx,dp),tsel(dx,dp),valsel(dx));
drhsc=arrayfun(@(i){@(x,p,dx,dp)drhs(i,x,p,dx,dp)},1:options.maxorder);
tau =    @(it,x,p)      dist_delay_tau(ind.msh.tau(it),bdsel(x,p),ind.msh,0);
dtau=@(ord,it,x,p,dx,dp)dist_delay_tau(ind.msh.tau(it),bdsel(dx,dp),ind.msh,ord);
dtauc=arrayfun(@(ord){@(it,x,p,dx,dp)dtau(ord,it,x,p,dx,dp)},1:options.maxorder);
funcs=set_funcs('wrap_rhs',rhs,...
    'wrap_dirderi',drhsc,...
    'wrap_tau',tau,...
    'wrap_dirdtau',dtauc,...
    'sys_ntau',@()ind.msh.ntau-1,...
    'sys_tau_seq',{ind.msh.tau},...
    'x_vectorized',true,'p_vectorized',true,...
    'xpattern',xpattern,...
    'lhs_matrix',zeros(ind.args.ival.nval,length(ixmap)));
%% append new funcs to tab
tab=dde_add_funcs(tab,options.name,funcs,...
    'par_old',ipmap,'x_old',ixmap,...
    'rhstau_args',ind.msh.ntau-1,...
    'nf',ind.args.ival.nval,pass_on{:});
end
%%
function [args,iunique]=local_from_global(args,name)
%% collect global indices
names=fieldnames(args);
nvars=length(names);
[indglob,types]=deal(cell(1,nvars)); %nvx stores array of variables type name
for i=1:nvars
    arg=args.(names{i});
    glob=[arg.global{:}];
    types{i}=find(strcmp([arg.type{:}],name));
    indglob{i}=glob(types{i});
end
%% map global to local (chain_delay specific) indices
indglobvec=[indglob{:}];
iunique=unique(indglobvec);
indloc=1:length(iunique);
loc_from_glob=NaN(1,max(iunique));
loc_from_glob(iunique)=indloc;
for i=1:nvars
    arg=args.(names{i});
    alocal=cell(1,length(arg.global));
    for k=1:length(alocal)
        istype=strcmp(arg.type{k},name);
        alocal{k}=NaN(1,length(arg.type{k}));
        alocal{k}(istype)=loc_from_glob([arg.global{k}(istype)]);
    end
    arg.local.(name)=alocal;
    arg.nval=length([arg.global{:}]);
    args.(names{i})=arg;
end
end
%%
function tau=dist_delay_tau(itau,bd,msh,order)
nvec=size(bd,2);
bd=reshape(bd,1,nvec);
if order<=1
    tau=reshape(msh.t(itau+1),[],1)*bd;
else
    tau=zeros(length(itau),nvec);
end
end
%%
function args = arg_wrap2cell(args)
argnames=fieldnames(args);
for i=1:length(argnames)
    arg=args.(argnames{i});
    if ~isempty(arg.replace)
        for k=1:length(arg.type)
            isrep=strcmp(arg.type{k},arg.replace{1});
            arg.type{k}(isrep)=repmat(arg.replace(2),1,sum(isrep));
        end
    end
    args.(argnames{i})=arg;
end
end
%%
function init=init_xpar_reduce(init)
init.local.parameter={[init.local.parameter{:}]};
init.local.x=        {[init.local.x{:}]};
init.type={[init.type{:}]};
init.global={[init.global{:}]};
end
%%
function ind_x=ind_local_x(arg)
nhist=length(arg.type);
ind_xc=cell(1,nhist);
for i=1:nhist
    ind_xc{i}=[arg.local.x{i}];
    pos_x=  strcmp(arg.type{i},'x');
    ind_xc{i}=ind_xc{i}(pos_x);
    ind_xc{i}=ind_xc{i}(:);
end
ind_x=cat(1,ind_xc{:});
end
%%
function xpval=int_delay_xparchoice(x,p,arg,ihist)
if nargin<4
    ihist=1;
end
ind_x=[arg.local.x{ihist}];
ind_par=[arg.local.parameter{ihist}];
nxp=length(arg.global{ihist});
isfun=true(1,nxp);
pos_x=  strcmp(arg.type{ihist},'x');
pos_par=strcmp(arg.type{ihist},'parameter');
isfun(pos_x|pos_par)=false;
nvec=size(p,3);
xpval=NaN(nxp,nvec);
ind_x=ind_x(pos_x);
ind_par=ind_par(pos_par);
xpval(pos_par,:)=reshape(p(1,ind_par,:),length(ind_par),nvec);
xpval(pos_x  ,:)=reshape(x(ind_x,  1,:),length(ind_x),  nvec);
ifun=find(isfun,1,'first');
if isempty(ifun) || ischar(arg.type{ihist}{ifun})
    return
end
xpval(isfun ,:)=arg.type{ihist}{ifun}(sum(isfun),nvec);
end
%%
function cfs=prod_deriv_coeffs(fac,mxord)
nfac=length(fac);
[cf,ct]=deal(cell(mxord,1));
cf{1}=ones(1,nfac);
ct{1}=1;
for k=2:mxord+1
    A=zeros(k*ones(1,nfac));
    szA=size(A);
    for i=1:size(cf{k-1},1)
        for j=1:nfac
            iv=cf{k-1}(i,:);
            iv(j)=iv(j)+1;
            ivc=num2cell(iv);
            ind=sub2ind(szA,ivc{:});
            A(ind)=A(ind)+1;
        end
    end
    sel=find(A(:)>0);
    cf{k}=NaN(length(sel),nfac);
    ct{k}=NaN(length(sel),1);
    for i=1:length(sel)
        iv=cell(1,nfac);
        [iv{:}]=ind2sub(szA,sel(i));
        cf{k}(i,:)=cat(2,iv{:});
        ct{k}(i)=A(sel(i));
    end
    tbrem=any(cf{k}>fac+1,2);
    cf{k}(tbrem,:)=[];
    ct{k}(tbrem)=[];
end
cfs.coeffs=cf;
cfs.count=ct;
end

function M=loc_set_matrix(M,ind)
nt=ind.args.t.nval;
nval=length(ind.args.ival.local.x{1});
if isnumeric(M) && isnan(M)
    M=speye(nt*length(ind.sum.input));
end
if iscell(M)
    M=reshape(cat(1,M{:}),nval,[]);
end
if isnumeric(M) && ndims(M)==3
    M=reshape(permute(M,[1,3,2]),nval,[]);
end
assert(size(M,2)==nt*length(ind.sum.input));
assert(size(M,1)==nval);
end