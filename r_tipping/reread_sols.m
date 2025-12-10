%% collected repeated calls when re-reading a solution from disk
function [prob,uidx,u0,maps]=reread_sols(prob,seglist,runlist,lablist,ip,varargin)
default={'match_plus_gamma',[],'add_gap_monitor',[],'fix_r_T',[],...
    'identify_parameters',[],'adjoint',''};
options=sco_set_options(default,varargin);
single_run=false;
if ~iscell(runlist)
    runlist=repmat({runlist},1,length(seglist));
    lablist=repmat({lablist},1,length(seglist));
    single_run=true;
end
adj=struct('do',~ismember(options.adjoint,{'','none'}),...
    'type',options.adjoint,'opt',[],'idx',[],'run',runlist{1},'lab',lablist{1}(end));
assert(~adj.do||single_run,...
    'adjoint computations can only be set up from single run, label');
for i=1:length(seglist)
    seg=seglist{i};
    bvp_args={seg,runlist{i},lablist{i}(end)};
    prob=ode_bvp2bvp(prob,bvp_args{:});
    [data.(seg),uidx.(seg),u0.(seg)]=...
        coco_get_func_data(prob,coco_get_id(seg,'bvp.seg1.coll'),'data','uidx','u0');
    maps.(seg)=data.(seg).pr.coll_seg.maps;
    ind=uidx.(seg);
    prob=coco_add_pars(prob,['T',seg],ind(maps.(seg).T_idx),...
        ['T',seg]);
    if ~adj.do
        continue
    end
    switch adj.type
        case 'init'
            prob=adjt_isol2bvp(prob,seg);
        case {'cont','switch'}
            prob=adjt_bvp2bvp(prob,bvp_args{:});
    end
    % this part should only apply for adjoint
    [adata, adj.idx.(seg)] = coco_get_adjt_data(prob,...
        coco_get_id(seg,'bvp.seg1.coll'), 'data', 'axidx');
    adj.opt.(seg) = adata.coll_opt;
    idx=adj.idx.(seg);
    opt=adj.opt.(seg);
    prob=loc_add_adjt(prob,['T',seg],adj,['d.T',seg],'aidx',idx(opt.T_idx));
end
if ischar(options.match_plus_gamma)
    prob=match_plus_gamma(prob,options.match_plus_gamma,uidx,maps,adj);
end
if ischar(options.identify_parameters)
    prob=identify_parameters(prob,options.identify_parameters,uidx,maps,ip,adj);
end
if ischar(options.add_gap_monitor)
    prob=add_gap_monitor(prob,options.add_gap_monitor,uidx,maps,adj);
end
if ischar(options.fix_r_T)
    prob=fix_r_T(prob,options.fix_r_T,uidx,u0,maps,ip,adj);
end
end
