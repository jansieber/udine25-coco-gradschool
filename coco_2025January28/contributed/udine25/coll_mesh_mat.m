%% Evaluate (derivative of) piecewise collocation polynomial
% assuming all requested points are inside the collocation interval
% mesh (assertion is placed)
%%
function [J,ix,tbp]=coll_mesh_mat(msh,t,varargin)
%% INPUT:
%
% * mshfine: 1xntstx(deg+1) sequence of times, where each ntst times
% correspond to subintervals on which solution is polynomial of degree deg
% if 1st arg is structure, we assume that it is coll_seg with entries
% mesh.tbp
% * ntst: number of subintervals, determines degree as numel(mshfine)/ntst-1
% * t (1 x nt): point(s) where to evaluate
%
%% Optional
%
% * loc of size(1,nt)): same size as t: if t(k)==msh(i) for some i in 2:N,
% then t fits to the interval msh([i-1,i]) if loc(k)==-1,
% otherwise, if loc(k)==+1, t fits to the interpolation in msh([i,i+1])
% For all other t loc is irrelevant. Default for loc is that for any i such
% that t(i)==t(i+1), loc(i) is -1, for all others it is 1.
% * kron (1): if kron(J,speye(n)) should be done at the end
% * sparse (true): matrix returned should be sparse
% * diff (integer, default 0): differentiate diff times
%
%% OUTPUT: 
%   J: sparse Jacobian:  P(y(:))(:)=J*y(:) for the piecewise collocation
%   polynomial on fullmesh(:)
%% default for loc
default={'kron',1,'diff',0,'sparse',true,'ix',[],'tbptol',1e-14};
[options,pass_on]=sco_set_options(default,varargin,'pass_on');
tbp = coll2msh(msh, options.tbptol);
[deg1,ntst]=size(tbp);
if numel(options.ix)<numel(t)
    [ix,t]=coll_mesh_find(tbp,t,pass_on{:});
else
    ix=options.ix;
end
nx=numel(t);
%% options
%% evaluate Lagrange polynomials on all interpolation times
dt=diff(tbp([1,end],:),[],1);
tscal=2*(t-tbp(1,ix))./dt(ix)-1;
%% evaluate barycentric weights
% for flexibility, (not assuming any particular interpolation grid)
submesh=2*(tbp(:,1)-tbp(1,1))/(tbp(end,1)-tbp(1,1))-1;
submesh([1,end])=[-1,1];
[w,Dmat]=barywt(submesh(:,1),options.diff);
%% calculate interpolation matrix
% using barycentric interpolation
ti_m=(ix-1)*(deg1)+1;
og=ones(deg1,1);
ox=ones(nx,1);
denom=tscal(og,:)-submesh(:,ox);
jac_ind(:,:,2)=ti_m(og,:)+repmat((0:deg1-1)',1,nx);
jac_ind(:,:,1)=repmat(1:nx,deg1,1);
fac=w(:,ox)./denom;
jac_vals=zeros(deg1,nx);
jac_vals(~denom(:))=1;
denomfin=all(denom~=0,1);
jac_vals(:,denomfin)=fac(:,denomfin)./repmat(sum(fac(:,denomfin),1),deg1,1);
%% Differentiate if requested
if options.diff>0
    div=dt(ix).^options.diff;
    jac_vals=(Dmat'*jac_vals)./repmat(div,deg1,1);
end
%% assemble Jacobian
jac_ind=reshape(jac_ind,[],2);
jac_vals=jac_vals(:);
J=sparse(jac_ind(:,1),jac_ind(:,2),jac_vals,nx,deg1*ntst);
if options.kron>1
    J=kron(J,speye(options.kron));
end
if ~options.sparse
    J=full(J);
end
end
%% construct barycentric weights for mesh
function [w,Dout]=barywt(submesh,difforder)
npi=length(submesh);
base_h=submesh';
base_v=submesh;
xdiff=base_h(ones(npi,1),:)-base_v(:,ones(npi,1));
w=1./prod(eye(npi)-xdiff,2);
Dout=eye(npi);
if difforder>0
    wrep=w(:,ones(npi,1));
    denom=(xdiff'-xdiff);
    denom(1:npi+1:end)=Inf;
    D=wrep'./wrep./denom;
    Drsum=sum(D,2);
    D(1:npi+1:end)=-Drsum;
    for i=1:difforder
        Dout=4*D*Dout;
    end
end
end