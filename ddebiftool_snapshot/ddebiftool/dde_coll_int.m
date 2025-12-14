function [y,J]=dde_coll_int(varargin)
%% integral of piecewise polynomial over a subinterval of mesh([1,end])
%% INPUT:
%
% * profile: profile on mesh t
% * mesh: representation points in [t1,t2] size(m*l+1)
% * x (2 x 1): interval  for integral: x(1) to x(2)
% * degree: order polynomials
%
% or
%
% * point: coll or psol structure
% * x (2 x 1): interval  for integral: x(1) to x(2)
%
%% OUTPUT: 
%	y integrals size(profile,1) x nx 
%   J: sparse Jacobian:  y(:)=J*profile(:)

%% distinguish between two calling modes numerical arguments
if isnumeric(varargin{1})
    assert(length(varargin)>=4);
    [profile,mesh,bd,degree]=deal(varargin{1:4});
    optind=5;
elseif isstruct(varargin{1})
    assert(length(varargin)>=2);
    [pt,bd]=deal(varargin{[1,2]});
    assert(isfield(pt,'profile'));
    [profile,mesh,degree]=deal(pt.profile,pt.mesh,pt.degree);
    optind=3;
else
    error('COLL:ARGS','dde_coll_eva: argument sequence not recognized');
end
%% options
default={'kron',false,'sparse',true,'output','profile','assert_boundaries',true};
options=dde_set_options(default,varargin(optind:end),'pass_on');
if options.assert_boundaries
    assert(all(bd(:)>=mesh(1)&bd(:)<=mesh(end)));
end
isneg=bd(1)>bd(2);
sgn=1-2*double(isneg);
if isneg
    bd=bd([2,1]);
end
y=[];
mesh_wt = dde_coll_int_mwt(mesh,degree,bd(:).');
int_mat=sgn*mesh_wt.wmat*mesh_wt.b2c;
if options.kron
    nx=size(profile,1);
    J=kron(int_mat,speye(nx));
else
    J=int_mat;
end
if strcmp(options.output,'profile') || nargout==2
    y=profile*int_mat.';
end
if strcmp(options.output,'matrix')
    [y,J]=deal(J,y);
end
end
