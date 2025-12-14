function E=dde_psol_combine(point,S,shift)
% only tested for shift with 0<=shift(:,1)<shift(:,2)
ns=size(S,1);
nx=size(S,2);
nd=size(shift,1);
t0=point.mesh(:);
nt=size(point.mesh,2);
E=sparse(ns*nt,nx*nt);
for i=1:nd
    ts=t0+shift(i,1)/shift(i,2);
    wrap=zeros(size(ts));
    wrap(ts<0)=floor(ts(ts<0));
    wrap(ts>1)=floor(ts(ts>1));
    ts(ts<0)=ts(ts<0)-wrap(ts<0);
    ts(ts>1)=ts(ts>1)-wrap(ts>1);
    E1=dde_coll_eva(point,ts(:)','output','matrix','kron',false);
    % the below is for rotations (where periodic bc is not true)
    E1(wrap~=0,end)=E1(wrap~=0,end)+wrap(wrap~=0);
    E1(wrap~=0,  1)=E1(wrap~=0,  1)-wrap(wrap~=0);
    E=E+kron(speye(nt),S(:,:,i))*kron(E1,speye(nx));
end
end