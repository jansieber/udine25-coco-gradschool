function funcs=dde_combined_funcs(tab,varargin)
%% define a function structure combined from function structures listed in tab
%% check for delayed derivatives
delayed_derivs=max(arrayfun(@(ref)ref.funcs.delayed_derivs,tab.fun_ref));
%% define combined rhs
sys_rhs=@(x,p)loc_rhs_dir(tab,delayed_derivs,0,x,p);
%% define combined drhs_dir, dtau_dir
drhs_dir=@(varargin)loc_rhs_dir(tab,delayed_derivs,varargin{:});
drhs_mf=@(varargin)loc_rhs_mf(tab,drhs_dir,varargin{:});
rhs_args={'wrap_rhs',sys_rhs,'drhs_dir',drhs_dir,'drhs_mf',drhs_mf};
%% check if we have variable delay anywhere
tp_delfcn=arrayfun(@(ref)ref.funcs.tp_del,tab.fun_ref);
tp_del=any(tp_delfcn);
if ~tp_del
    itaupar=loc_sys_tau_const(tab);
    sys_tau=@()itaupar;
    assert(tab.ntau==length(itaupar));
    sys_tau_seq={1:tab.ntau};
    tau_args={'sys_tau',sys_tau,'sys_tau_seq',sys_tau_seq};
else
    dtau_dir=@(varargin)loc_dtau(tab,varargin{:});
    dtau_mf=@(varargin)mult_deriv(dtau_dir,{1,[2,2]},varargin{:},'nf',tab.ntau,'splitcomplex',true);
    fseqlist=arrayfun(@(ref){ref.funcs.sys_tau_seq},tab.fun_ref);
    sys_tau_seq=join_sys_tau_seq(tab.fun_ref,fseqlist);
    sys_tau=@(it,x,p)loc_sys_tau_tpdel(it,x,p,tab);
    tau_args={'wrap_tau',sys_tau,...
        'dtau_dir',dtau_dir,'dtau_mf',dtau_mf,...
        'sys_ntau',@()tab.ntau,'sys_tau_seq',sys_tau_seq};
end
%% define combined lhs_matrix
lhs_matrix=loc_lhsmatrix(tab);
%% find how many derivatives were provided
dirderi_num=min(arrayfun(@(f)f.funcs.dirderi_provided(),tab.fun_ref));
%% combine xpatterns
xpattern=loc_xpattern(tab);
%% all functions have been wrapped already
funcs=set_funcs('lhs_matrix',lhs_matrix,...
    rhs_args{:},...
    tau_args{:},...
    'delayed_derivs',delayed_derivs,...
    'dirderi_num',dirderi_num,...
    'xpattern',xpattern,...
    varargin{:});
end
%%
function y=loc_rhs_mf(tab,drhs_dir,varargin)
incl_deriv=~isempty(varargin) && ischar(varargin{1}) &&...
    strcmp(varargin{1},'incl_deriv');
dims={1,[2+double(incl_deriv),2]};
y=mult_deriv(drhs_dir,dims,varargin{1+double(incl_deriv):end},'nf',tab.nf,'splitcomplex',true);
end
%%
function y=loc_rhs_dir(tab,delayed_derivs,varargin)
incl_deriv=~isempty(varargin) && ischar(varargin{1}) &&...
    strcmp(varargin{1},'incl_deriv');
iarg=double(incl_deriv);
deriv_needed=delayed_derivs+1;
assert(incl_deriv||deriv_needed==1);
[ord,x,p,iarg]=deal(varargin{iarg+(1:3)},iarg+3);
if incl_deriv
    [nx,ntau,nderiv]=size(x,[1,2,3]);
else
    [nx,ntau]=size(x,[1,2]);
    nderiv=1;
end
assert(deriv_needed<=nderiv);
x=reshape(x,nx,ntau,nderiv,[]);
vecdim=[size(x,4:ndims(x)),1];
nvec=size(x,4);
p=reshape(p,1,size(p,2),[]);
y=NaN(tab.nf,nvec);
for i=1:length(tab.fun_ref)
    fcn=tab.fun_ref(i);
    ixtau=[1,1+fcn.rhstau_args];
    [nxi,nti,ndxi]=deal(length(fcn.x),length(ixtau),fcn.funcs.delayed_derivs+1);
    xshape=[{nxi,nti},repmat({ndxi},1,double(ndxi>1)),{[]}];
    xi=reshape(x(fcn.x,ixtau,1:ndxi,:),xshape{:});
    pari=p(1,fcn.par,:);
    if ord==0
        y(fcn.f,:)=fcn.funcs.wrap_rhs(xi,pari);
    else
        dx=reshape(varargin{iarg+1},nx,ntau,nderiv,[]);
        dxi=reshape(dx(fcn.x,ixtau,1:ndxi,:),xshape{:});
        dpari=varargin{iarg+2}(1,fcn.par,:);
        y(fcn.f,:)=fcn.funcs.drhs_dir(ord,xi,pari,dxi,dpari);
    end
end
y=reshape(y,[size(y,1),vecdim]);
end
%%
function y=loc_dtau(tab,ord,x,p,dx,dp)
[nx,ntau]=size(x,[1,2]);
if nargin<=4
    [dx,dp]=deal(NaN(size(x)),NaN(size(p)));
end
[x,dx]=deal(reshape(x,nx,ntau,[]),reshape(dx,nx,ntau,[]));
vecdim=[size(x,3:ndims(x)),1];
nvec=size(x,3);
p=reshape(p,1,size(p,2),[]);
y=NaN(ntau-1,nvec);
for i=1:length(tab.fun_ref)
    fcn=tab.fun_ref(i);
    fun=fcn.funcs.dtau_dir;
    ixtau=[1,1+fcn.rhstau_args];
    xi=x(fcn.x,ixtau,:);
    pari=p(1,fcn.par,:);
    dxi=dx(fcn.x,ixtau,:);
    dpari=dp(1,fcn.par,:);
    tau_inc0=fun(ord,xi,pari,dxi,dpari);
    y(fcn.tauprovided,:)=tau_inc0(2:end,:);
end
y=[zeros(1,nvec);y];
y=reshape(y,[size(y,1),vecdim]);
end
%%
function y=loc_sys_tau_tpdel(it,x,p,tab)
[nx,ntau]=size(x,[1,2]);
vecdim=[size(x,3:ndims(x)),1];
x=reshape(x,nx,ntau,[]);
nf=length(it);
if nf==0
    y=zeros([0,vecdim]);
    return
end
f_ind=unique(tab.taumap(it));
y=NaN(nf,size(x,3));
for i=1:length(f_ind)
    ref=tab.fun_ref(f_ind(i));
    [iftau,itauprov]=intersect(it,ref.tauprovided);
    [dum,ind]=ismember(iftau,ref.rhstau_args); %#ok<ASGLU>
    ixtau_all=[1,1+ref.tauprovided];
    ixtau=ixtau_all(1:min(ind));
    xarg=x(ref.x,ixtau,:);
    pararg=p(1,ref.par,:);
    y(itauprov,:)=ref.funcs.wrap_tau(ind,xarg,pararg);
end
y=reshape(y,[size(y,1),vecdim]);
end
%%
function lhs_matrix=loc_lhsmatrix(tab)
lhs_matrix=zeros(tab.nf,length(tab.xtest));
nfun=length(tab.fun_ref);
for i=1:nfun
    ref=tab.fun_ref(i);
    lhs_matrix(ref.f,ref.x)=ref.funcs.lhs_matrixfun(length(ref.x));
end
end
%%
function itau=loc_sys_tau_const(tab)
nfun=length(tab.fun_ref);
itauc=arrayfun(@(ref){reshape(ref.funcs.sys_tau(),1,[])},tab.fun_ref);
itau=NaN(1,tab.ntau);
for i=1:nfun
    itau(tab.fun_ref(i).tau)=tab.fun_ref(i).par(itauc{i});
end
assert(all(~isnan(itau)));
end
%% define combined sys_tau_seq
function jointseq=join_sys_tau_seq(fun_ref,fseqlist)
nfun=length(fun_ref);
for i=1:nfun
    if isempty(fun_ref(i).tauprovided)
        fseqlist{i}={};
    else
        fseqlist{i}=cellfun(@(k){reshape(fun_ref(i).tauprovided(k),1,[])},...
            fseqlist{i});
    end
end
seqlen=cellfun(@length,fseqlist);
if all(seqlen==0)
    jointseq=cell(1,0);
    return
end
fseqlist=fseqlist(seqlen>0);
tauord=cellfun(@(s){cat(1,reshape([s{:}],1,[]),dde_tauin_from_tauseq(s))},fseqlist);
tmin=cellfun(@(t)t(1,1),tauord);
[dum,is]=sort(tmin); %#ok<ASGLU>
tauord=tauord(is);
tauord_comb=reshape(cat(2,tauord{:}),2,[]);
tauord_comb(2,:)=tauord_comb(2,:)+...
    cumsum([0,-min(diff(tauord_comb(2,:)),0)]);
jointseq=dde_tauin_from_tauseq(tauord_comb(2,:),'reverse');
end
%% define combined xpattern
function xpattern=loc_xpattern(tab)
xpattern=[];
for i=1:length(tab.fun_ref)
    ref=tab.fun_ref(i);
    ix=ref.x;
    itau=[1,1+ref.rhstau_args];
    xpat_i=ref.funcs.xpattern;
    xpat_i(1,:)=ix(xpat_i(1,:));
    xpat_i(2,:)=itau(xpat_i(2,:));
    xpattern=dde_join_xpattern(xpat_i,xpattern);
end
end
