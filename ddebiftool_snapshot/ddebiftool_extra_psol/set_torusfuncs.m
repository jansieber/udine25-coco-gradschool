function [trfuncs,extra_freepar,initfuncs]=set_torusfuncs(funcs,point,method,varargin)
%% set up extended systems, numbering and values of additional parameters and artificial delays
%% default extra conditions
standardcond={'psol_phase_condition',true,'Torus_norm',true,'Torus_phase_condition',true};
%% process options
default={'biftype','torus','nullparind',zeros(0,1),...
    'usercond',cell(0,1),'initcond',cell(0,1),standardcond{:}}; %#ok<CCAT>
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% set up numbering and values of additional parameters
ip.dim=size(point.profile,1);    % dimension of original problem
ip.extdim=2*ip.dim;
ip.xrg=1:ip.dim;
ip.re=ip.dim+(1:ip.dim);
ip.im=ip.re(end)+(1:ip.dim);
ip.null=[ip.re,ip.im];
ip.nuserpar=length(point.parameter); % number of original system parameters
ip.omega=ip.nuserpar+1;             % location of add. parameter omega
ip.period=ip.omega+1;            % location of add. parameter (equal to period)
ip=dde_ip_addnullpar(ip,point,options.nullparind,0,'lastname','omega');
%% set up functions of extended system
trfuncs=set_psolbiffuncs(funcs,ip,@sys_rhs_TorusBif_var,method,pass_on{:});
trfuncs.kind=options.biftype;
trfuncs.get_comp=@(p,component)extract_from_tr(p,component,options.biftype,ip);
torus_temp=dde_torus_template(trfuncs,point);
%% additional free parameters
extra_freepar=[ip.omega,ip.period];
%% add requested and standard extra conditions
% for initialization and for continuation
usercond=dde_test_cond(options.usercond,torus_temp);
initcond=dde_test_cond(options.initcond,torus_temp);
embeddedconds=dde_embed_cond(trfuncs,funcs.sys_cond,'solution');
trfuncs=dde_funcs_add_cond(trfuncs,embeddedconds,'name','Torus_embedded');
initfuncs=dde_funcs_add_cond(trfuncs,initcond,'name','Torus_initcond');
trfuncs=dde_funcs_add_cond(trfuncs,usercond,'name','Torus_usercond');
psol_phase_cond=dde_sys_cond_create('name','psol_phasecondition',...
    'fun',@(p,pref)sys_cond_psol_phase_condition(p,pref,ip.dim),'reference',true);
extra_norm_cond=dde_sys_cond_create('name','Torus_norm',...
    'fun',@(p)sys_cond_Torus_norm(p,ip.dim),'reference',false);
extra_phase_cond=dde_sys_cond_create('name','Torus_phase_condition',...
    'fun',@(p,pref)sys_cond_Torus_phase_condition(p,pref,ip.dim),'reference',true);
sys_cond_extra=cat(1,psol_phase_cond,extra_phase_cond,extra_norm_cond);
fixperiod=sys_cond_coll_fixperiod('fixperiod',ip.period);
sys_extra_sel=sys_cond_extra(...
    [options.psol_phase_condition;options.Torus_phase_condition;options.Torus_norm]);
trfuncs=dde_funcs_add_cond(trfuncs,[fixperiod;sys_extra_sel]);
sys_init_sel=sys_cond_extra(options.psol_phase_condition);
initfuncs=dde_funcs_add_cond(initfuncs,sys_init_sel);
end