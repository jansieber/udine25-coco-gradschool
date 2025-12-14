function [branch,succ,fail,rjct]=br_contn(funcs,branch,max_tries,varargin)
%% extend DDE-BIFTOOL branch
% function [c_branch,succ,fail,rjct]=br_contn(funcs,branch,max_tries)
% INPUT:
%   * funcs problem functions
%	* branch initial branch (contains method and previous points)
%	* max_tries maximum number of tries
%   * optional (named):
%     'plotaxis' (default []): if plotting is on
%     (branch.method.continuation.plot>=1) then the user may specify an axis
%     into which to plot, default [] chooses gca
%   * other optional arguments are passed on as fields of method subfields
%   
% OUTPUT:
%	branch extended branch
%	succ number of succesfull corrections
%	fail number of failed corrections
%	rjct number of failures which ended in rejecting a point
%
% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
%% introduce optional argument to set plotting axis
branch=replace_branch_pars(branch,branch.parameter.free,varargin);
branch=br_bound_delay(funcs,branch); %% bound delays
method=branch.method;
free_par=branch.parameter.free;
%% install check for upper and lower bounds for parameters 
% (mthcont is modified method.continuation with standard stops added
mthcont=br_install_standardfcns(funcs,branch,varargin{:});
%% Initialize counters and logical variables
[stat.fail,stat.rjct,predstop,corstop]=deal(0);
stat.max_tries=max_tries;
udata=[];
old_success=1;
assert(length(branch.point)>1,...
    'br_contn:start','BR_CONTN: could not start branch, length=%d.',length(branch.point));
l_usetangent=mthcont.use_tangent && isfield(branch.point(end),'nvec');
emptypoint=repmat(branch.point(1),1,0);
steplength=p_norm(p_axpy(-1,branch.point(end-1),branch.point(end)));
%% loop performing maximally max_tries steps
for tries=1:max_tries
    stat.tries=tries;
    assert(length(branch.point)>1,'br_contn:fail',['BR_CONTN: could not continue branch,',...
        '%d points, %d fails, %d rejected.'],tries-stat.fail,stat.fail,stat.rjct);
    %% find secant, add tangent if required (condition checked inside add_tangent)
    secant=p_axpy(-1,branch.point(end-1),branch.point(end));
    extrapoint=emptypoint;
    if ~old_success && ~l_usetangent % if previous step failed ,track back
        extrapoint=branch.point(end);
        branch.point=branch.point(1:end-1);
    end
    branch.point(end)=add_tangent(funcs,method,branch.point(end),secant,free_par,varargin{:});
    last_point=branch.point(end);
    %% predict and determine steplength
    dist=p_norm(secant);
    factor=(1/2)*double(~old_success)+mthcont.steplength_growth_factor*double(old_success);
    if l_usetangent
        pred_dir=last_point.nvec.tangent.vector;
        scprod=p_dot(secant,pred_dir);
        br_warn(mthcont,{'warn_angle',scprod<mthcont.warn_angle*dist},...
            'br_contn: angle between tangent and secant<%3.0fdeg\n',acosd(mthcont.warn_angle));
        steplength=factor*steplength; % shrink distance if angle sharp
    else
        pred_dir=p_axpy(1/dist,secant,[]);
        steplength=factor*dist;
    end
    %% break if steplength falls below minimum (if set)
    if isfield(mthcont,'steplength_minimum') && abs(steplength)<mthcont.steplength_minimum
        br_warn(mthcont,{},'br_contn: stepsize falls below minimum %g\n',...
            mthcont.steplength_minimum);
        break
    end
    pred_point=p_axpy(steplength,pred_dir,last_point);
    %% check for maximal steplengths
    [pred_point,steplength]=check_maxstep(pred_point,steplength,...
        branch.parameter.max_step,last_point,pred_dir);
    %% choose reference point for correction r.h.s
    ref_point=last_point;
    %% check if user has instructed to stop, plot or print
    test_point=pred_point;
    [predstop,udata]=br_online_userfcn('predictor',udata,mthcont,branch,stat,test_point);
    if predstop>0
        break
    end
    %% correct: choose pseudoarclength condition
    stpcond=p_secant(pred_dir,p_norm(ref_point));
    if ~mthcont.steplength_condition
        stpcond=[];
    end
    [cor_point,cor_success]=...
        p_correc(funcs,pred_point,free_par,stpcond,method.point,tries+1,ref_point,varargin{:});    
    pdiff=p_axpy(-1,ref_point,cor_point);
    %% check if outcome is finite
    ndiff=p_norm(pdiff);
    cor_success=cor_success && isfinite(ndiff);
    %% compute angle between points
    % if minimal angle requested: call too large angle unsuccessful
    % take into account that after previously not successful step the new
    % point is between last and previous point
    cospdiff=p_dot(pdiff,pred_dir);
    angle_bad=l_usetangent && isfield(mthcont,'minimal_angle') && ...
        cospdiff<mthcont.minimal_angle*ndiff;
    br_warn(mthcont,{'warn_angle',~angle_bad,cospdiff<mthcont.warn_angle*ndiff},...
        'br_contn point %d: angle between predictor and new secant small\n',length(branch.point)+1);
    cor_success=cor_success && ~angle_bad;
    %% accumulate failures and removed points, put cor_success and angle into stat
    stat.fail=stat.fail+double(~cor_success);
    stat.rjct=stat.rjct+double(~cor_success && ~old_success && ~l_usetangent);
    [stat.success,stat.cosang,stat.steplength]=deal(cor_success,cospdiff/ndiff,steplength);
    if cor_success
        %% check if user has instructed to stop, plot or print
        test_point=cor_point;
        [corstop,udata]=br_online_userfcn('corrector',udata,mthcont,branch,stat,test_point);
        if corstop>0
            break
        end
        %% ensure format fits, insert point & previous last point (if present)
        cor_point=dde_trim_point(cor_point,branch.point(1));
        append=[cor_point,extrapoint];
        branch.point(end+(1:length(append)))=append;
    elseif ~old_success && mthcont.halt_before_reject~=0
            break
    end
    old_success=cor_success;
end
%% final corrections and actions after continuation
[branch,stat]=br_final_userfcn(corstop,predstop,mthcont,udata,branch,stat,[ref_point,test_point]);
succ=stat.tries-stat.fail;
[rjct,fail]=deal(stat.rjct,stat.fail);
end
%% install stop functions
function mth=br_install_standardfcns(funcs,branch,varargin)
method=branch.method;
mth=branch.method.continuation;
free_par=branch.parameter.free;
%% add stop checks
mth=br_add_stop(mth); % convert method.continuation.stop field to struct of correct format
%% install check for upper and lower bounds for parameters
parcheck=@(pt,output)dde_check_parameterbounds(funcs,method,pt,free_par,...
    branch.parameter.min_bound,branch.parameter.max_bound,output);
mth=br_add_stop(mth,'name','parameterbounds',...
    'online',@(p)parcheck(p,'flag'),'final',@(p)parcheck(p,'point'),'state','predictor');
mth=br_add_stop(mth,'name','parameterbounds',...
    'online',@(p)parcheck(p,'flag'),'final',@(p)parcheck(p,'point'),'state','corrector');
%% for variable delay problems, install check if delay hits zero
if funcs.tp_del
    delcheck=@(pt,output)dde_check_delaysign(funcs,method,pt,free_par,output);
    mth=br_add_stop(mth,'name','delayzero',...
        'online',@(p)delcheck(p,'flag'),'final',@(p)delcheck(p,'point'),'state','predictor');
end
%% install online plotting if required
if branch.method.continuation.plot_progress>0
    mth=br_add_stop(mth,'name','onlineplot_pred','argument_format','branch','state','predictor',...
        'online',@(data,branch,stats,point)br_online_plot('predictor',data,branch,stats,point,varargin{:}));
    mth=br_add_stop(mth,'name','onlineplot_corr','argument_format','branch','state','corrector',...
        'online',@(data,branch,stats,point)br_online_plot('corrector',data,branch,stats,point,varargin{:}));
end
%% install online printing if required
if isfield(mth,'print_progress') && mth.print_progress
    mth=br_add_stop(mth,'name','onlineprint','argument_format','branch','state','corrector',...
        'online',@(data,branch,stats,point)br_online_print('corrector',data,branch,stats));
end
end
%% check for maximal steplengths
function [new_point,steplength]=check_maxstep(new_point,steplength,max_step,last_point,pred_dir)
fraction=1;
for j=1:size(max_step,1)
    if max_step(j,1)==0
        dp=p_norm(p_axpy(-1,last_point,new_point));
    elseif max_step(j,1)==-1
        if isfield(new_point,'period')
            lp=setfield(last_point,'period',new_point.period); %#ok<SFLD>
        else
            lp=last_point;
        end
        dp=p_norm(p_axpy(-1,lp,new_point));
    else
        dp=abs(new_point.parameter(max_step(j,1))-last_point.parameter(max_step(j,1)));
    end
    if dp>max_step(j,2)
        f=max_step(j,2)/dp;
        if f<fraction
            fraction=f;
        end
    end
end
if fraction<1
    steplength=steplength*fraction;
    new_point=p_axpy(steplength,pred_dir,last_point);
end
end
%% add tangent
function point=add_tangent(funcs,mth,point,secant,free_par,varargin)
if ~mth.continuation.use_tangent || ~isfield(point,'nvec')
    return
end
if ~isfield(point.nvec,'tangent') || ...
        ~length(point.nvec.tangent.free_par)==length(free_par) || ...
        ~all(point.nvec.tangent.free_par==free_par)
    stpcond=p_tangent(funcs,mth.point,point,free_par,'border',secant,varargin{:});
    stpcond=p_axpy(1/p_norm(stpcond),stpcond,[]);
    point.nvec.tangent=struct('free_par',free_par,'vector',stpcond);
end
cosangpdiff=p_dot(point.nvec.tangent.vector,secant);
if cosangpdiff<0
    point.nvec.tangent.vector=p_axpy(-1,point.nvec.tangent.vector,[]);
end
end
