%% check sign of delays
function varargout=dde_check_delaysign(funcs,mth,pts,free_par,output)
dostop=false;
crit_delay=struct('nr',[],'t_z',[]);
if funcs.tp_del~=0 && ...
        (~isfield(mth.continuation,'permit_negative_delay') || ...
        ~mth.continuation.permit_negative_delay)
    [crit_delay.nr,crit_delay.t_z]=p_tsgn(funcs,pts(end));
    if ~isempty(crit_delay.nr)
        dostop=true;
        if isfield(mth.continuation,'warnings') && ...
                mth.continuation.warnings && strcmp(output,'flag')
            fprintf('delay %d crossed zero\n',crit_delay.nr);
        end
    end
end
if strcmp(output,'flag')
    varargout={dostop};
    return
end
if ~dostop || length(pts)<2
    varargout={pts(end),false};
    return
end
%% construct initial guess
tau_n=p_tau(funcs,pts(end),crit_delay.nr,crit_delay.t_z);
tau_p=p_tau(funcs,pts(end-1),crit_delay.nr,crit_delay.t_z);
del_tau=tau_p-tau_n;
pp=p_axpy(-tau_n/del_tau,pts(end-1),[]);
pts(end)=p_axpy(tau_p/del_tau,pts(end),pp);
%% correct final point
pini=pts(end);
npar=length(pini.parameter);
itz=npar+(1:length(crit_delay.t_z));
pini.parameter(itz)=crit_delay.t_z;
free_par=[free_par,itz];
del_cond=@(pt,pref)dde_apply({'dde_',pini.kind,'_delay_zero_cond'},pt,funcs,...
    free_par,crit_delay.nr,itz);
fcn=dde_funcs_add_cond(funcs,del_cond,'name','delaysign','reference',true);
mth.point.extra_condition=true;
pref=pts(end-1);
pref.parameter(itz)=crit_delay.t_z;
[pout,suc]=p_correc(fcn,pini,free_par,[],mth.point,0,pref);
pout.parameter(itz)=[];
varargout={pout,suc};
end
%%