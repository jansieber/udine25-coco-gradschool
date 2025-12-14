function x_d=dde_select_tau(x,taudim,derivdim,rhs_xderiv)
%% return x(1:nx,1:ntau,1:nderivs+1) from x w format n1 x n2 x n3 x n4 x ...
sz=size(x);
maxdim=max([length(sz),taudim,derivdim]);
sz=[sz,ones(1,maxdim-length(sz))];
nderivs=max(rhs_xderiv)+1;
otherdims=setdiff(1:length(sz),[taudim,derivdim]);
x=reshape(permute(x,[otherdims,taudim,derivdim]),[],sz(derivdim));
x_d=reshape(x(:,1:nderivs),[sz(otherdims),sz(taudim)*nderivs]);
perm=1:length(sz)-1;
perm([length(perm),taudim])=[taudim,length(perm)];
x_d=permute(x_d,perm);
end
