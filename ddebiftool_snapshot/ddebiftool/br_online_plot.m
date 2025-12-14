function [flag,data,branch]=br_online_plot(state,data,branch,stats,new_point,varargin)
%% Online plotting during continuation
% this function is called during |br_contn|.
%
% *Inputs*
%
% * |state|: flag ='corrector' or 'predictor', indicating where we are
% * |ax|: axis object where plot is updated
% * |branch|: branch structure maintained during continuation. Field
% |'method.continuation'| contains fields relevant for plotting such as
% |plot_measure|.
% * |new_point|: predicted or corrected point (end point of line thatis
% drawn)
% * |last_point|: previous/reference point (starting point of line)
% * |step|: integer, counting the step (used in title for user information)
% * |props|: properties, eg, color used for line (hard coded in |br_contn| at the moment),
% different colors indicate predictor step and correced step.
% * |txt|: replacement text used for printout if
% method.continuation.plot_progress>0 but <1 (hard coded as |'pred'| and |'
% cor'| in |br_contn| for now).
%%
flag=false;
if branch.method.continuation.plot_progress==0
    return
end
if isempty(data) && branch.method.continuation.plot_progress>=1
    data.ax=gca;
end
default={{'plotaxis','ax'},data.ax,...
    'predictor_prop',{'color',[0.85, 0.325, 0.098],'linewidth',0.5},...
    'corrector_prop',{'color',[0,    0.447, 0.741],'linewidth',1,'MarkerSize',9},...
    'predictor_txt','pred',...
    'corrector_txt','             cor',...
    'point_prop',{'color',[0,0,0],'linewidth',2,'MarkerSize',5,'Marker','+'}};
options=dde_set_options(default,varargin,'pass_on');
data.ax=options.plotaxis;
if ~any(strcmp(state,{'predictor','corrector'}))
    return
end
props=options.([state,'_prop']);
txt=options.([state,'_txt']);
last_point=branch.point(end);
if length(new_point)>1
    new_point=new_point(end);
end
step=length(branch.point);
mth=branch.method.continuation;
if mth.plot<=0
    return
end
m_isfun=false;
if isempty(mth.plot_measure)
    [x_m,y_m]=df_measr(0,branch);
elseif iscell(mth.plot_measure)
    x_m=mth.plot_measure{1};
    y_m=mth.plot_measure{2};
    m_isfun=true;
else
    x_m=mth.plot_measure.x;
    y_m=mth.plot_measure.y;
end
if ~isempty(x_m) && isstruct(x_m)
    x1=p_measur(last_point,x_m);
    x2=p_measur(new_point,x_m);
elseif m_isfun
    x1=feval(x_m,last_point);
    x2=feval(x_m,new_point);
else
    x1=step;
    x2=step+1;
end
if ~isempty(y_m) && isstruct(y_m)
    y1=p_measur(last_point,y_m);
    y2=p_measur(new_point,y_m);
elseif m_isfun
    y1=feval(y_m,last_point);
    y2=feval(y_m,new_point);
else
    y1=step;
    y2=step+1;
end
if mth.plot>=1
    try
        ish=ishold(data.ax);
        hold(data.ax,'on');
        if ~dde_isoctave()
            args={'HandleVisibility','off'};
        else
            args={};
        end
        plot(data.ax,[x1 x2],[y1 y2],props{:},args{:});
        plot(data.ax,x2,y2,'.',props{:},args{:});
        if strcmp(state,'corrector')
            if isfield(data,'pthandle') && ishghandle(data.pthandle)
                set(data.pthandle,'XData',x2,'YData',y2);
            else
                data.pthandle=plot(data.ax,x2,y2,options.point_prop{:});
            end
        end
        steplength=p_norm(p_axpy(-1,last_point,new_point));
        title(data.ax,sprintf('%d points, steplen=%3.1e, fail:%d, rjct:%d',...
            step,abs(steplength),stats.fail,stats.rjct));
        if ~ish
            hold(data.ax,'off');
        end
    catch ME
        warning('br_contn:online_plot',...
            'br_contn: online plotting requested but plot command failed.\n%s',ME.message);
    end
    if mth.plot_progress
        drawnow;
    end
elseif mth.plot_progress
    fprintf('%s: (%g, %g)\n',txt,x2,y2);
end
end
