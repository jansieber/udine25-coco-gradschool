function [pfuncs,extra_freepar,initfuncs]=set_POfoldfuncs(funcs,point,method,varargin)
%% set up extended systems, numbering and values of additional parameters and artificial delays
%% default extra conditions
standardcond={'psol_phase_condition',true,'POEV1_norm',true,'POEV1_phase_condition',true};
%% process options
default={'nullparind',zeros(0,1),...
    'usercond',cell(0,1),'initcond',cell(0,1),standardcond{:}}; %#ok<CCAT>
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
ip.dim=size(point.profile,1);
ip.extdim=ip.dim;
ip.nuserpar=length(point.parameter); % number of original system parameters
%% extend problem
ip.beta=ip.nuserpar+1;       % location of add. parameter beta
[ip,extrapar]=dde_ip_addnullpar(ip,point,options.nullparind,1,...
    'lastname','beta','periodnullpar',ip.beta);
pfuncs=set_psolbiffuncs(funcs,ip,@sys_rhs_POEV1_var,method,pass_on{:});
pfuncs.usermethod=method;
pfuncs.kind='POfold';
get_comp=@(p,component)extract_from_POEV1(p,component,ip);
pfuncs.get_comp=get_comp;
pfold_temp=dde_pofold_template(pfuncs,point);
%% additional free parameters
extra_freepar=[ip.beta,extrapar];
%% add requested and standard extra conditions
% for initialization and for continuation
usercond=dde_test_cond(options.usercond,pfold_temp);
initcond=dde_test_cond(options.initcond,pfold_temp);
embeddedconds=dde_embed_cond(pfuncs,funcs.sys_cond,'solution');
pfuncs=dde_funcs_add_cond(pfuncs,embeddedconds,'name','POfold_embedded');
initfuncs=dde_funcs_add_cond(pfuncs,initcond,'name','POfold_initcond');
pfuncs=dde_funcs_add_cond(pfuncs,usercond,'name','POfold_usercond');
psol_phase_condition=dde_sys_cond_create('name','psol_phasecondition',...
    'fun',@(p,pref)sys_cond_psol_phase_condition(p,pref,ip.dim),'reference',true);
POEV1_phase_condition=dde_sys_cond_create('name','POEV1_phase_condition',...
    'fun', @(p,pref)sys_cond_POEV1_phase_condition(p,pref,ip.dim),'reference',true);
POEV1_norm=dde_sys_cond_create('name','POEV1_norm','fun',...
    @(p)sys_cond_POEV1_norm(p,ip,extra_freepar,'res',1,'period',false),...
        'reference',false);
sys_cond_extra=cat(1,psol_phase_condition,POEV1_phase_condition,POEV1_norm);
sys_extra_sel=sys_cond_extra(...
    [options.psol_phase_condition;options.POEV1_phase_condition;options.POEV1_norm]);
pfuncs=dde_funcs_add_cond(pfuncs,sys_extra_sel);
initfuncs=dde_funcs_add_cond(initfuncs,sys_extra_sel);
end
