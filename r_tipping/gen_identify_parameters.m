%% equate parameters from two bvp or coll segments
function prob=gen_identify_parameters(prob,gluename,uidx,maps,varargin)
default={'free',[]};
opts=sco_set_options(default,varargin,'pass_on');
npars=length(maps(1).p_idx);
include=setdiff(1:npars,opts.free);
for i=2:-1:1
    gu_idx{i}=uidx{i}(maps(i).p_idx(include));
end
prob=coco_add_glue(prob,gluename,gu_idx{:});
end
