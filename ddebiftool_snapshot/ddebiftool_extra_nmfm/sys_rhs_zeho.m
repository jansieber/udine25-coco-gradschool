function res=sys_rhs_zeho(funcs,order,xx,p,dxx,dp)
%% r.h.s. and derivative for Zero-Hopf bifurcation
%
%%
n=size(xx,1)/4;
sq=@(a)reshape(a,n,[]);
ntau=dde_num_delays(funcs)+1;
ntaurep=@(x)repmat(reshape(x,n,1,[]),1,ntau,1);
[ix,iq0,iq1r,iq1i]=deal(1:n,n+(1:n),2*n+(1:n),3*n+(1:n));
[ipu,iom]=deal(1:size(p,2)-1,size(p,2));
[xu,pu,q0,q1r,q1i]=deal(sq(xx(ix,1,:)),p(1,ipu,:),...
    sq(xx(iq0,1,:)),sq(xx(iq1r,1,:)),sq(xx(iq1i,1,:)));
q1=q1r+1i*q1i;
omega=reshape(p(1,iom,:),1,[]);
Delta0=@(varargin)ch_matrix(funcs,xu,pu,0,       varargin{:});
Delta1=@(varargin)ch_matrix(funcs,xu,pu,1i*omega,varargin{:});
rhs=@(order,varargin)funcs.drhs_dir(order,ntaurep(xu),pu,varargin{:});
switch order
    case 0
        ruser=rhs(0);
        rfold=Delta0('dx',q0);
        rhopf=Delta1('dx',q1);
    case 1
        [dxu,dpu,dq0,dq1r,dq1i]=deal(sq(dxx(ix,1,:)),dp(1,ipu,:),...
            sq(dxx(iq0,1,:)),sq(dxx(iq1r,1,:)),sq(dxx(iq1i,1,:)));
        dq1=dq1r+1i*dq1i;
        dom=repmat(reshape(dp(1,iom,:),1,[]),n,1);
        ch_dxp_args={'devx',dxu,'devp',dpu,'order',1};
        ruser=rhs(1,ntaurep(dxu),dpu);
        rfold=Delta0('dx',dq0)+Delta0('dx',q0, ch_dxp_args{:});
        rhopf=Delta1('dx',dq1)+Delta1('dx',q1,'deri',1)*1i.*dom+...
              Delta1('dx', q1, ch_dxp_args{:});
end
res=cat(1,ruser,rfold,real(rhopf),imag(rhopf));
end
