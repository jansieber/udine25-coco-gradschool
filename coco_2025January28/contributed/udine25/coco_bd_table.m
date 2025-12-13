function bd=coco_bd_table(run,varargin)
default={'cut',0,'numlab',false};
opts=sco_set_options(default,varargin,'pass_on');
bdc=coco_bd_read(run);
names=bdc(1,1:end-opts.cut);
vals=bdc(2:end,1:end-opts.cut);
[~,ix]=unique(names);
is=sort(ix);
bd=cell2table(vals(:,is),'VariableNames',names(is));
if ~opts.numlab
    return
end
labs=bd{:,'LAB'};
lab_has_val=cellfun(@(c)~isempty(c),labs);
numlabs=NaN(length(labs),1);
numlabs(lab_has_val)=cat(1,labs{lab_has_val});
bd.LAB=numlabs;
end
