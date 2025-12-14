function [extfuncs,extra_freepar,initfuncs]=set_foldfuncs(funcs,point,varargin)
%% set up extended systems, numbering and values of additional parameters
%% process options
default={'hjac',1e-4,'nullparind',zeros(0,1),...
    'usercond',cell(0,1),'initcond',cell(0,1)};
options=dde_set_options(default,varargin,'pass_on');
ip.dim=size(point.x,1);
ip.extdim=ip.dim;
ip.nuserpar=length(point.parameter); % number of original system parameters
%% extend problem
[ip,extra_freepar]=dde_ip_addnullpar(ip,point,options.nullparind,1);
extfuncs=funcs;
extfuncs.kind='fold';
get_comp=@(p,component)extract_from_fold(p,component,ip);
extfuncs.ip=ip;
extfuncs.get_comp=get_comp;
par0=point.parameter;
par0(extra_freepar)=1;
fold_temp=dde_fold_create('point',point,'parameter',par0);
%% add requested and standard extra conditions
% for initialization and for continuation
usercond=dde_test_cond(options.usercond,fold_temp);
initcond=dde_test_cond(options.initcond,fold_temp);
extfuncs.sys_cond=repmat(dde_sys_cond_create(),0,1);
embeddedconds=dde_embed_cond(extfuncs,funcs.sys_cond,'solution');
embed=@(p,component,template)dde_stst_extendblanks(p,template);
extfuncs.embed=embed;
extfuncs=dde_funcs_add_cond(extfuncs,embeddedconds,'name','fold_embedded');
initfuncs=dde_funcs_add_cond(extfuncs,initcond,'name','fold_initcond');
extfuncs=dde_funcs_add_cond(extfuncs,usercond,'name','fold_usercond');
end
%% extract components from fold solution branch or point array
function result_array=extract_from_fold(fold_array,comp,ip)
%
%
%% check if input is branch rather than point array
if ~isfield(fold_array,'kind') && isfield(fold_array,'point')
    fold_array=fold_array.point;
end
%% extract named components
dim=ip.dim;
npar=ip.nuserpar;
type={'kind','solution','nullvector','solution_for_stability',...
    'eigenvector'};
for i=1:length(fold_array)
    fold=fold_array(i);
    switch comp
        case type{1} %'kind'
            result_array='fold';
            break
        case type{2} %'solution'
            result=fold;
            result.x=result.x(1:dim,:);
            result.parameter=result.parameter(1:npar);
        case type([3,5]) %'nullvector','eigenvector'
            result=fold.v;
            result.parameter=zeros(1,npar);
            result.parameter(ip.nullparind(:,1))=result.parameter(ip.nullparind(:,2));
        case type{4} % including flag
            result=fold;
            result.x=result.x(1:dim,:);
            result.v=result.v(1:dim,:);
            result.parameter=result.parameter(1:npar);
            result.flag='fold';
        otherwise
            fprintf('known component types:\n');
            for k=1:length(type)
                fprintf('%s\n',type{k});
            end
            result_array=[];
            break;
    end
    result_array(i)=result; %#ok<AGROW>
end
end
