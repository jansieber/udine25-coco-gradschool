function mesh_wt = dde_coll_int_mwt(tmesh,deg,bd,varargin)
%% create weights and nodes for integration (Gauss-Legendre)
% from bd(1) to bd(2), bd(2) to bd(3) ... to bd(end) for increasing bd
% on each subinterval of joint mesh. mesh_wt.t_c
% is sequence of nodes, mesh_wt.w is weight for each node, mesh_wt.b2c is
% map from mesh nodes to collocation nodes t_c, mesh_wt.b2bd is map from
% mesh nodes to boundary bd(end:-1:1)
default={'t_quad',[],'w_quad',[],'mesh_wt',[]};
options=dde_set_options(default,varargin,'pass_on');
if ~isempty(options.mesh_wt)
    mesh_wt=options.mesh_wt;
    return
end
if isempty(options.t_quad)
    [t_quad,w_quad]=legpts(deg,[0,1]);
else
    [t_quad,w_quad]=deal(options.t_quad,options.w_quad);
end
nbd=length(bd);
tcoarse=tmesh(1:deg:end);
[tcoarse,~,ibd]=unique([bd,tcoarse]);
tcoarse=tcoarse(ibd(1):ibd(nbd));
w=dde_coll_meshfill(tcoarse,1,'grid',w_quad,'acc',false);
nw=length(w);
t_c=dde_coll_meshfill(tcoarse,1,'grid',t_quad);
b2c=dde_coll_eva([],tmesh,t_c,deg,'kron',false,'output','matrix');
b2bd=dde_coll_eva([],tmesh,bd,deg,'kron',false,'output','matrix');
irow=interp1(bd(:),[1:nbd-1,nbd-1]',t_c(:),'previous')';
wmat=sparse(irow(:),(1:nw)',w,nbd-1,nw);
mesh_wt=struct('t_c',t_c,'w',w,'wmat',wmat,'b2c',b2c,'b2bd',b2bd);
end