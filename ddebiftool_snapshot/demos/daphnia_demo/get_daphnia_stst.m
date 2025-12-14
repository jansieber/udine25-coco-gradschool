function u0=get_daphnia_stst(par0,ix,ip,pname,pval,r_rg)
if nargin<6
    r_rg=[0,2];
end
par0(ip.(pname))=pval;
ststfun=dde_sym2fun(@sym_ststeq,'rhs');
ststcomp=dde_sym2fun(@sym_stst,'rhs');
nx=length(fieldnames(ix));
u0=NaN(nx,1);
u0(ix.r)=fzero(@(x)ststfun(x,par0)-1,r_rg);
u0(2:end)=ststcomp(u0(ix.r),par0);
end
