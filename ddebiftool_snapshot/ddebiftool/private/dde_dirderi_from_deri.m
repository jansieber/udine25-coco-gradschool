function J=dde_dirderi_from_deri(x,p,dx,dp,funcs,itau)
%% wrapper to create dirderi from old style deri type derivatives
xsize=size(x);
psize=size(p);
vecdim=[xsize(3:end),1];
x=reshape(x,xsize(1),xsize(2),[]);
p=reshape(p,psize(1),psize(2),[]);
dx=reshape(dx,size(x));
dp=reshape(dp,psize(1),psize(2),[]);
if nargin<=6 || isempty(itau)
    deri=funcs.wrap_deri_old;
    ntau=xsize(2);
    nf=size(funcs.lhs_matrixfun(xsize(1)),1);
else
    deri=@(xa,pa,nxa,npa,v)funcs.wrap_dtau_old(itau,xa,pa,nxa,npa);
    ntau=min([xsize(2),itau]);
    nf=length(itau);
end
nvec=size(x,3);
J=zeros(nf,nvec);
xsel=find(any(any(dx,1),3));
xsel=xsel(xsel<=ntau);
psel=find(any(any(dp,1),3));
%xsel=1:ntau;
%psel=1:psize(2);
for j=xsel
    Avec=reshape(deri(x(:,1:ntau,:),p,j-1,[],[]),nf,[]);
    dxvec=sparse_blkdiag(dx(:,j,:));
    A_d=Avec*dxvec;
    J=J+A_d;
end
for j=psel
    Avec=reshape(deri(x(:,1:ntau,:),p,[],j,[]),nf,[]);
    dpvec=sparse_blkdiag(dp(1,j,:));
    A_d=Avec*dpvec;
    J=J+A_d;
end
J=reshape(J,[nf,vecdim]);
end
