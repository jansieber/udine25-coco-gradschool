function Delta=dde_stst_ch_matrix(A,tau,lambda,varargin)
%% combine matrices Aj to characteristic matrix Delta(lambda)
% and  its derivatives
%
%%
default={'deri',0,'dxp',false,'lhs_matrix',eye(size(A,1))};
options=dde_set_options(default,varargin,'pass_on');
lfac=[lambda,1,zeros(1,options.deri-1)];
[nf,nx,r]=size(A,1:3);
max_rhs_xderiv=size(A,4)-1;
if length(tau)<r
    taus=[0;tau(:)];
else
    taus=tau(:);
end
taus=taus(:,ones(1,max_rhs_xderiv+1));
explt=exp(-lambda*taus);
for k=0:max_rhs_xderiv
    explt(:,k+1)=(-taus(:,k+1)).^(k+options.deri).*explt(:,k+1);
end
explt=repmat(reshape(explt,1,1,r,[]),nf,nx,1,1);
if ~options.dxp
    Delta_ini=options.lhs_matrix;
else
    Delta_ini=zeros(nf,nx);
end
Delta = lfac(options.deri+1)*Delta_ini-sum(sum(A.*explt,3),4);
end
