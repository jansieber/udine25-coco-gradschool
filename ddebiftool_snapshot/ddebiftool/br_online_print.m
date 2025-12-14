function [flag,data,branch]=br_online_print(state,data,branch,stats,varargin)
flag=false;
mth=branch.method.continuation;
if ~isfield(mth,'print_progress') || ~mth.print_progress||~strcmp(state,'corrector')
    return
end
if stats.success
    suc='success';
else
    suc='reject ';
end
angstr=sprintf(', cos(ang)=%g',stats.cosang);
fprintf('tries: %d of %d, pt %d %s steplength=%g%s\n',...
    stats.tries,stats.max_tries,length(branch.point)+1,suc,stats.steplength,angstr);
end

