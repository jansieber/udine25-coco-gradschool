%% collected repeated calls when re-reading a solution from disk
function [prob,uidx,u0,maps]=reread_sols(prob,seglist,runlist,lablist,iv,ip,varargin)
default={'match_plus_gamma',[],'add_gap_monitor',false,'fix_r_T',[],...
    'identify_parameters',[],'adjoint','','prepend',{}};
[options,pass_on]=sco_set_options(default,varargin,'pass_on');
prefix=options.prepend;
do_copy=~isempty(prefix);
nopars={};
if do_copy
    nopars={'-no-pars'};
end
single_run=false;
if ~iscell(runlist)
    runlist=repmat({runlist},1,length(seglist));
    lablist=repmat({lablist},1,length(seglist));
    single_run=true;
end
assert(~do_copy||single_run,...
    'restart can only be set up from single run, label');
for i=1:length(seglist)
    [seg,run,lab]=deal(seglist{i},runlist{i},lablist{i}(end));
    name=seg;
    segname=seg;
    if do_copy 
        if ~strncmp(options.prepend{1},seg,length(options.prepend{1}))
            name=coco_get_id(options.prepend{:},seg);
        else
            segname=seg(length(options.prepend{1})+2:end);
        end
    end
    bvp_args={name,run,seg,lab,nopars{:}}; %#ok<CCAT>
    prob=ode_bvp2bvp(prob,bvp_args{:});
    [data.(segname),uidx.(segname),u0.(segname)]=...
        coco_get_func_data(prob,coco_get_id(name,'bvp.seg1.coll'),'data','uidx','u0');
    maps.(segname)=data.(segname).pr.coll_seg.maps;
    ind=uidx.(segname);
    prob=coco_add_pars(prob,coco_get_id(prefix{:},['T',segname]),ind(maps.(segname).T_idx),...
            coco_get_id(prefix{:},['T',segname]));
end
if ischar(options.match_plus_gamma)
    prob=match_plus_gamma(prob,prepend(options,'match_plus_gamma'),iv,uidx,maps);
end
ids=options.identify_parameters;
for i=1:length(ids)
    names={ids(i).match1,ids(i).match2};
    gluename=coco_get_id(options.prepend{:},sprintf('%s_%s_glue',names{:}));
    prob=gen_identify_parameters(prob,gluename,{uidx.(names{1}),uidx.(names{2})},...
        [maps.(names{1}),maps.(names{2})],'free',ids(i).free);
end
if ischar(options.add_gap_monitor)
    prob=add_gap_monitor(prob,prepend(options,'add_gap_monitor'),iv,uidx,maps,pass_on{:});
end
if ischar(options.fix_r_T)
    prob=fix_r_T(prob,prepend(options,'fix_r_T'),uidx,u0,maps,ip);
end
end
%%
function newname=prepend(options,name)
newname=options.(name);
if isempty(newname)
    return
end
newname=coco_get_id(options.prepend{:},options.(name));
end