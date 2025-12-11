%% equate parameters from different segments to each other
function prob=identify_parameters(prob,gluename,uidx,maps,ip,varargin)
default={'free_names',{}};
opts=sco_set_options(default,varargin,'pass_on');
npars=length(maps.u_minus.p_idx);
include=1:npars;
for i=1:length(opts.free_names)
    include(ip.(opts.free_names{i}))=[];
end
u1_idx=[uidx.u_minus(maps.u_minus.p_idx);...
        uidx.u_gamma(maps.u_gamma.p_idx(include))];
u2_idx=[uidx.u_plus( maps.u_plus.p_idx);...
        uidx.u_plus( maps.u_plus.p_idx(include))];
prob=coco_add_glue(prob,gluename,u1_idx,u2_idx);
end
