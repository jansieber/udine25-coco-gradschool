function [info,atlas,bd]=info_from_run(run,varargin)
%% read name info for atlas (requires adding query func to problem)
atlas=coco_bd_read(run,'atlas');
bd=coco_bd_read(run,'bd');
bdinfo=coco_bd_read(run,'bddat');
info.pnames=bdinfo.bdp_names;
info.nbh=cellfun(@(c)c.neigh,atlas.charts,'uniformoutput',false);
%info.id=cellfun(@(c)c.id,atlas.charts);
%info.nbh=info.nbh(info.id);
info.p=cell2mat(cellfun(@(c)c.p,atlas.charts,'uniformoutput',false))';
%info.p=info.p(info.id,:);
info.xp=cell2mat(cellfun(@(c)c.xp,atlas.charts,'uniformoutput',false))';
%info.xp=info.xp(info.id,:);
[info.ip2x,info.x2pdist]=sco_match_vector(info.xp,info.p,varargin{:});
info.xpnames=info.pnames(info.ip2x);
info.npts=length(info.nbh);
info.dim=atlas.charts{1}.dim;
end
