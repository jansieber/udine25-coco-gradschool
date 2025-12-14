function [prob,uidx,u0,maps]=reread_sols_copy(prob,run,lab,iv,ip,varargin)
default={'h_dev',1e-3,'fix_r_T',false,'init',true};
opts=sco_set_options(default,varargin,'pass_on');
id_pars=@(name1,name2,free)struct('match1',name1,'match2',name2,'free',free);
glue_segments=[id_pars('um','up',[]),id_pars('ug','up',ip.phi)];
seglist={'ug','up','um'};
if opts.fix_r_T
    args={'fix_r_T','rtfix'};
else
    args={};
end
[prob,uidx{1},u0{1},maps(1)]=reread_sols(prob,seglist,run,lab(1),iv,ip,...
    'match_plus_gamma','pglue','add_gap_monitor','gap','identify_parameters',glue_segments,...
    args{:});
prep={'prepend',{'copy'}};
if opts.init
    seg2=seglist(1,:);
else
    seg2=cellfun(@(s){coco_get_id('copy',s)},seglist);
end
[prob,uidx{2},u0{2},maps(2)]=reread_sols(prob,seg2,run,lab(end),iv,ip,prep{:},...
    'match_plus_gamma','pglue','add_gap_monitor','gap','identify_parameters',glue_segments);
if opts.fix_r_T
    prob=coco_add_glue(prob,'glue_Tcopies',...
        [uidx{1}.up(maps(1).up.T_idx),uidx{1}.um(maps(1).um.T_idx)],...
        [uidx{2}.up(maps(2).up.T_idx),uidx{2}.um(maps(2).um.T_idx)]);
end
prob=gen_identify_parameters(prob,'glue_copies',{uidx{1}.up,uidx{2}.up},...
    [maps(1).up,maps(2).up],'free',[ip.r,ip.phi]);
prob=coco_add_glue(prob,'diff_phifix',uidx{1}.um(maps(1).um.p_idx(ip.phi)),...
    uidx{2}.um(maps(2).um.p_idx(ip.phi)),opts.h_dev);
prob=coco_add_functionals(prob,'diff_r',[1,-1]/opts.h_dev,0,[uidx{1}.um(maps(1).um.p_idx(ip.r));...
    uidx{2}.um(maps(2).um.p_idx(ip.r))],'dr','inactive');
end