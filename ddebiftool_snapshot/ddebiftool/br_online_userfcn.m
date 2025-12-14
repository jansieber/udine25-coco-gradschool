function [dostop,udata]=br_online_userfcn(state,udata,mth,branch,stats,point)
%% check user monitoring functions such as parameter bound checks delay=zero checks
% online plotting, printing or saving
dostop=0;
if ~isfield(mth,'stops')
    return
end
if isempty(udata)
    udata=repmat({[]},length(mth.stops),1);
end
pts=[branch.point(end),point];
argslist=struct(...
    'branch',@(data){data,branch,stats,pts},...
    'points',@(data){pts});
calls=struct(...
    'branch',@(f,ind,args)f(args{:}),...
    'points',@(f,ind,args)deal(f(args{:}),udata{ind}));
selection=find(arrayfun(@(s)strcmp(s.state,state),mth.stops));
for i=1:length(selection)
    ind=selection(i);
    form=argslist.(mth.stops(ind).argument_format);
    args=form(udata{ind});
    call=calls.(mth.stops(ind).argument_format);
    [lstop,udata{ind}]=call(mth.stops(ind).online,ind,args);
    if lstop
        if ~dostop
            dostop=ind;
        end
        br_warn(mth,{},'br_contn warning: boundary %d hit: %s.\n',...
            ind,mth.stops(dostop).name);
    end
end
end
