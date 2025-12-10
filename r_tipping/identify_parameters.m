%% equate parameters from different segments to each other
function prob=identify_parameters(prob,gluename,uidx,maps,ip,adj)
npars=length(maps.u_minus.p_idx);
nonphi=setdiff(1:npars,ip.phi);
u1_idx=[uidx.u_minus(maps.u_minus.p_idx);...
        uidx.u_gamma(maps.u_gamma.p_idx(nonphi))];
u2_idx=[uidx.u_plus( maps.u_plus.p_idx);...
        uidx.u_plus( maps.u_plus.p_idx(nonphi))];
prob=coco_add_glue(prob,gluename,u1_idx,u2_idx);
if nargin<6 ||~adj.do
    return
end
a1_idx=[adj.idx.u_minus(adj.opt.u_minus.p_idx);...
        adj.idx.u_gamma(adj.opt.u_gamma.p_idx(nonphi))];
a2_idx=[adj.idx.u_plus( adj.opt.u_plus.p_idx);...
        adj.idx.u_plus( adj.opt.u_plus.p_idx(nonphi))];
prob=loc_add_adjt(prob,gluename,adj,'aidx',[a1_idx(:);a2_idx(:)]);
end
