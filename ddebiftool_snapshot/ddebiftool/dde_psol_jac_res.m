%% residual & jacobian for psol
%
% $Id: dde_psol_jac_res.m 369 2019-08-27 00:07:02Z jansieber $
%%
function varargout=dde_psol_jac_res(funcs,pt,free_par,method,varargin)
default={{'collocation_parameters','c'},method.collocation_parameters,...
    'phase_condition',method.phase_condition,'rotationcheck',true,...
    'bcmat',eye(size(pt.profile,1)),'pref',[],'phase_tolerance',1e-8,...
    'output','J','append_tc',false};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
pt=dde_coll_check(pt);
if isfield(method,'maxsize')
    pass_on=[{'maxsize'},{method.maxsize},pass_on(:)'];
end
[Js,res]=dde_coll_jac(funcs,pt,free_par,...
    'c',options.collocation_parameters,'matrix',method.matrix,...
    'Dtmat',funcs.lhs_matrixfun(size(pt.profile,1)),...
    'output',options.output,pass_on{:}); % beware that for output=='res', Js=res
if strcmp(options.output,'res')
    [Js,res]=deal(res,Js);
end
neqs=length(res);
ieq.de=1:neqs;
[ind,len]=dde_ind_from_point(pt,free_par);
%% periodicity condition:
resbc=pt.profile(:,1)-pt.profile(:,end);
if options.rotationcheck
    resbc=mod(resbc+pi,2*pi)-pi;
end
resbc=resbc(:);
dresbc=dde_coll_bd_diff(pt,'c',options.collocation_parameters,...
    'boundary',true,pass_on{:});
resbc=options.bcmat*cat(2,resbc,dresbc);
[nbc,nx]=deal(numel(resbc),size(pt.profile,1));
res=[res;resbc(:)];
ieq.bc=neqs+(1:nbc);
Jbc(:,ind.profile(:,1))=eye(nx);
Jbc(:,ind.profile(:,end))=-eye(nx);
dJbc=dde_coll_bd_diff(pt,...
    'c',options.collocation_parameters,'boundary',true,'output','matrix',pass_on{:});
Jbcmat=reshape(options.bcmat*reshape(cat(1,Jbc,dJbc),nx,[]),[],numel(pt.profile));
%% phase condition:
[resph,Jph_dx]=loc_phasecond(pt,ind,options);
ieq.ph=ieq.bc(end)+(1:length(resph));
res=[res;resph];
%% append info about smoothness conditions
ieq.de_f=ieq.de(1:end-Js.nforce_smooth);
ieq.de_smooth=ieq.de(end-Js.nforce_smooth+1:end);
if options.append_tc
    ieq.tc=Js.tc;
end
%% output for res
varargout=dde_setupoutput('jac_res',options.output,sparse(length(res),len),res,ieq);
if strcmp(options.output,'res')
    return
end
%% assemble Jacobian
Js.profile=cat(1,Js.profile,Jbcmat);
Js.profile=[Js.profile;Jph_dx];
J(:,ind.profile)=sparse(Js.profile);
J(1:neqs,ind.period)=Js.period;
J(1:neqs,ind.parameter)=Js.parameter;
if strcmp(method.matrix,'full')
    J=full(J);
end
varargout=dde_setupoutput('jac_res',options.output,J,res,ieq);
end

function [resph,Jph_dx] = loc_phasecond(pt,ind,options)
if ischar(options.phase_condition) ||options.phase_condition
    if isempty(options.pref)
        pref=pt;
        resph_to_0=true;
    else
        pref=dde_coll_check(options.pref);
        resph_to_0=false;
        if max(abs(diff(pref.profile,[],2)))<options.phase_tolerance
            resph_to_0=true;
        end
    end
    [resph,Jph_pt]=sys_cond_psol_phase_condition(pt,pref);
    Jph_dx(1:ind.profile(end))=Jph_pt.profile(:).';
    if resph_to_0
        resph=0;
    end
else
    resph=zeros(0,1);
    Jph_dx=sparse(0,ind.profile(end));
end
end