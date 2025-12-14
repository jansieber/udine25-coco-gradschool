function [branch,stat]=br_final_userfcn(corstop,predstop,mth,udata,branch,stat,pts)
%% final corrections and actions after continuation
dostop=corstop+predstop;
if dostop>0 && ~isempty(mth.stops(dostop).final) && strcmp(mth.stops(dostop).argument_format,'points')
        %% continuation stopped at boundary and correction installed
        [newpoint,bd_success]=mth.stops(dostop).final(pts);
        if bd_success
            branch.point(end+1)=newpoint;
        else
            stat.fail=stat.fail+1;
        end
elseif corstop>0
    branch.point(end+1)=pts(end);
end
order=arrayfun(@(s)s.order,mth.stops);
force=arrayfun(@(s)s.force,mth.stops);
present=arrayfun(@(s)~isempty(s.final),mth.stops);
todo=find(force & present);
[dum,ix]=sort(order(todo)); %#ok<ASGLU>
todo=todo(ix);
for i=1:length(todo)
    [dum,dum,branch]=mth.stops(todo(i)).final(udata{todo(i)},branch,stat,pts); %#ok<ASGLU>
end
end
