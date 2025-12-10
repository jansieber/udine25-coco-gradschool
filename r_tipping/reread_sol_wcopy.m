%% collected repeated calls when re-reading a solution from disk
function [prob,uidx,u0,maps]=reread_sol_wcopy(prob,seglist,runlist,lablist,ip,varargin)
default={'match_plus_gamma',[],'add_gap_monitor',[],'fix_r_T',[],...
    'identify_parameters',[],'adjoint','','prepend',{}};
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
    [seg,run,lab]=deal(seglist{i},runlist{i},lablist{i}(end));
    bvp_args={seg,run,lab};
    if isempty(options.prepend)
        prob=ode_bvp2bvp(prob,bvp_args{:});
        [data.(seg),uidx.(seg),u0.(seg)]=...
            coco_get_func_data(prob,coco_get_id(seg,'bvp.seg1.coll'),'data','uidx','u0');
        maps.(seg)=data.(seg).pr.coll_seg.maps;
        ind=uidx.(seg);
        prob=coco_add_pars(prob,['T',seg],ind(maps.(seg).T_idx),...
            ['T',seg]);
    else
        prefix=options.prepend;
        name=coco_get_id(options.prepend{:},seg);
        [bvpsol,bvpdata.(seg)]=bvp_read_solution(bvp_args{:});
        [collsol,colldata]=coll_read_solution(coco_get_id(seg,'bvp.seg1'),'r_of_phi',2)
        frhs={colldata.fhan,colldata.dfdxhan,colldata.dfdphan};
        pnames=
        prob=ode_isol2bvp(prob,name,frhs{:},collsol.tbp,collsol.xbp,collsol.p,umnames,bc_um);
        
    end
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
    prob=match_plus_gamma(prob,prepend(options,'match_plus_gamma'),uidx,maps,adj);
end
if ischar(options.identify_parameters)
    prob=identify_parameters(prob,prepend(options,'identify_parameters'),uidx,maps,ip,adj);
end
if ischar(options.add_gap_monitor)
    prob=add_gap_monitor(prob,prepend(options,'add_gap_monitor'),uidx,maps,adj);
end
if ischar(options.fix_r_T)
    prob=fix_r_T(prob,prepend(options,'fix_r_T'),uidx,u0,maps,ip,adj);
end
end
%%
function newname=prepend(options,name)
newname=coco_get_id(options.prepend{:},options.(name));
end