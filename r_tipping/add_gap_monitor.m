%% introduce equation and variable monitoring gap between x_-(1) and x_+(0)
function prob=add_gap_monitor(prob,fname,uidx,maps,varargin)
default={'indicator','GAP','gap_parname',fname};
opts=sco_set_options(default,varargin,'pass_on');
if ~isempty(opts.gap_parname)
    gap_pars={opts.gap_parname,'inactive'};
else
    gap_pars={};
end
prob=coco_add_glue(prob,fname,...
    uidx.u_plus(maps.u_plus.x0_idx(1)),...
    uidx.u_minus(maps.u_minus.x1_idx(1)),...
    0,gap_pars{:});
if ~isempty(opts.gap_parname)
    prob=coco_add_event(prob,opts.indicator,opts.gap_parname,0);
end
end
