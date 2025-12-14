%% equate parameters from different segments to each other
function prob=identify_parameters(prob,gluename,uidx,maps,ip,varargin)
default={'free_names',{}};
opts=sco_set_options(default,varargin,'pass_on');
npars=length(maps.um.p_idx);
include=1:npars;
for i=1:length(opts.free_names)
    include(ip.(opts.free_names{i}))=[];
end
u1_idx=[uidx.um(maps.um.p_idx);...
        uidx.ug(maps.ug.p_idx(include))];
u2_idx=[uidx.up( maps.up.p_idx);...
        uidx.up( maps.up.p_idx(include))];
prob=coco_add_glue(prob,gluename,u1_idx,u2_idx);
end
