function branch=br_add_save(branch,name,varargin)
default={'folder',[pwd(),filesep,'data',filesep],'batchsize',1,'cleanup',true,...
    'partname','_pt','onlineargs',{},'finalargs',{}};
options=dde_set_options(default,varargin,'pass_on');
if ~(exist(options.folder,'dir')==7)
    mkdir(options.folder);
end
options.name=name;
options.basept=[options.folder,filesep,name,options.partname];
options.ext='.mat';
options.baseptname=@(srg)[options.basept,srg,options.ext];
options.basefull=[options.folder,filesep,name,options.ext];
options.onlinelist={};
branch=br_add_stop(branch,'name','save','state','predictor','argument_format','branch',...
    'online',@(data,branch,stats,point)br_online_save('predictor',data,branch,stats,options),...
    'final',@(data,branch,stats,point)br_final_save(data,branch,options.finalargs),...
    'force',true,'order',2);
end
function [flag,data]=br_online_save(state,data,branch,stats,options)
flag=false;
if ~strcmp(state,'predictor')
    return
end
if isempty(data)
    data=options;
    data.numlen=floor(1+log10(length(branch.point)+stats.max_tries));
    data.formatstr=['%0',num2str(data.numlen),'d'];
end
if mod(length(branch.point),data.batchsize)==0
    ind=length(branch.point)-data.batchsize+1;
    point=branch.point(ind:end);
    fname=data.baseptname(num2str(ind,data.formatstr));
    save(fname,'point',options.onlineargs{:});
    data.onlinelist{end+1}=fname;
end
end
function [flag,data,branch]=br_final_save(data,branch,finalargs)
flag=false;
s=struct(data.name,branch);
save(data.basefull,'-struct','s',finalargs{:});
if data.cleanup
    delete(data.onlinelist{:});
end
end
