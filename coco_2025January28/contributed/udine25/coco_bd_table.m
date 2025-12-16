function bd=coco_bd_table(run,varargin)
default={'cut',0,'numlab',false};
opts=sco_set_options(default,varargin,'pass_on');
if ~iscell(run)
    bdc=coco_bd_read(run);
else
    bdc=run;
end
names=bdc(1,1:end-opts.cut);
vals=bdc(2:end,1:end-opts.cut);
[~,ix]=unique(names);
is=sort(ix);
try
    bd=cell2table(vals(:,is),'VariableNames',names(is));
catch % old matlab version: table headers need to be valid variable names
    varnames=matlab.lang.makeValidName(names(is));
    bd=cell2table(vals(:,is),'VariableNames',varnames);
    bd.Properties.VariableDescriptions=names(is);
end    
if ~opts.numlab
    return
end
labs=bd{:,'LAB'};
lab_has_val=cellfun(@(c)~isempty(c),labs);
numlabs=NaN(length(labs),1);
numlabs(lab_has_val)=cat(1,labs{lab_has_val});
bd.LAB=numlabs;
end
