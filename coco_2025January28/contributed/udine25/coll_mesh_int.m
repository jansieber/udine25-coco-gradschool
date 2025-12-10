function [inth,L,U,P,Q,S,J]=coll_mesh_int(msh,x)
msh=coll2msh(msh);
[degp1,ntst]=size(msh);
[nt,deg]=deal(degp1*ntst,degp1-1);
lgbase=legpts(deg,[0,1]);
tlg=lgbase(:)*diff(msh([1,end],:),[],1)+msh(ones(deg+1,1),:);
J=[coll_mesh_mat(msh,[1;tlg(:)]');sparse(ntst-1,nt)];
D=coll_mesh_mat(msh,tlg(:)','diff',1);
ont=ones(1,ntst);
C=sparse(1:ntst,1:degp1:nt,ont,ntst,nt)+...
    sparse(2:ntst,degp1:deg+1:nt-1,-ont(2:end),ntst,nt);
Jint_inv=cat(1,C(1,:),D,C(2:end,:));
[L,U,P,Q,S]=lu(Jint_inv);
if nargin>1
    intproj=coll_mesh_mat(msh,x(:).');
    %intproj=reshape(permute(cat(3,-proj(1:2:end,:),proj(2:2:end,:)),[1,3,2]),[],numel(x));
    inth=intproj*Q*(U\(L\(P*(S\J))));
else
    inth=@(x)x*Q*(U\(L\(P*(S\J))));
end
end