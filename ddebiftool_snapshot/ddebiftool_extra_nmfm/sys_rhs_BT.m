function res=sys_rhs_BT(funcs,x_ext,p,dx_ext,dp)
%% r.h.s. for Takens-Bogdanov bifurcation
% sys_cond_BT needs to be appended
%
% Delta(0)v=-dxf([x,x,..],p)[v,v...]
% [d/dlam]Delta(0)v=M*v+dxf([x,x,...],p)[0*v,tau1*v,...]
%%
[nx,nvec]=deal(size(x_ext,1)/3,size(x_ext,3));
sq=@(x)reshape(x,nx,[]);
[ir,irq1,irq2]=deal(1:nx,nx+(1:nx),2*nx+(1:nx));
[iv1,iv2]=deal(1:nvec,nvec+(1:nvec));
[x,q1,q2]=deal(sq(x_ext(ir,:,:)),sq(x_ext(irq1,:,:)),sq(x_ext(irq2,:,:)));
ntau=dde_num_delays(funcs)+1;
ntaurep=@(x)repmat(reshape(x,nx,1,[]),1,ntau,1);
Delta=@(varargin)ch_matrix(funcs,x,p,0,varargin{:});
Dq =Delta('dx',[q1,q2]);
dDq=Delta('dx', q1,   'deri',1);
if nargin<4 %r.h.s.
    res=cat(1,funcs.wrap_rhs(ntaurep(x),p),...
          Dq(:,iv1),...
         dDq(:,iv1)+Dq(:,iv2));
    return
end
%% 1st derivative
[dx,dq1,dq2]=deal(sq(dx_ext(ir,:,:)),sq(dx_ext(irq1,:,:)),sq(dx_ext(irq2,:,:)));
Ddq   =Delta('dx',[dq1,dq2]);
dDdq  =Delta('dx', dq1,     'deri',1);
Dq_xp =Delta('dx',[ q1, q2],           'devx',[dx,dx],'devp',cat(3,dp,dp),'order',1);
dDq_xp=Delta('dx',  q1,     'deri',1,  'devx', dx,    'devp',      dp,    'order',1);
res=cat(1,funcs.drhs_dir(1,ntaurep(x),p,ntaurep(dx),dp),...
     Dq_xp(:,iv1)+ Ddq(:,iv1),...
    dDq_xp(:,iv1)+dDdq(:,iv1)+Dq_xp(:,iv2)+Ddq(:,iv2));
end
