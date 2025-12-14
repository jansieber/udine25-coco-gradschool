function varargout=p_correc_rhs(funcs,method,point,free_par,varargin)
%% residual and jacobian of nonlinear equation to be solved by p_correc
% if optional argument x is given then point stucture is only used as
% template. Otherwise, it is used to extract x using dde_x_from_point.
% Optional point structures in step_cond are added as rows to the Jacobian.
% Some nonlinearities require a reference point, given as optional pref.
%
% $Id: p_correc_rhs.m 348 2019-06-19 13:09:01Z jansieber $
%%
default={'x',[],...
    'step_cond',repmat(point,0,1),...
    'pref',repmat(point,0,1),...
    'extracolumns',[],...
    'output','res+J'};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% extract extra parameters if needed
extracolumns=options.extracolumns;
[extralen,nextrapar]=size(extracolumns);
%% if given use optional arg x as variable
if ~isempty(options.x)
    point=dde_point_from_x(options.x(1:end-nextrapar),point,free_par);
    extrapar=options.x(end-nextrapar+1:end);
else
    extrapar=point.parameter(free_par(end-nextrapar+1:end));
    free_par=free_par(1:end-nextrapar);
end
%% add extra conditions as required
res_extra=[];
cond_extra=repmat(point,0,1);
if method.extra_condition
    [res_extra,cond_extra,ieq.extra_names]=call_sys_cond(funcs,point,options.pref);
end
%% get residual and jacobian of r.h.s
out=cell(1,3);
[out{:}]=dde_apply({'dde_',point.kind,'_jac_res'},...
    funcs,point,free_par,method,'pref',options.pref,pass_on{:},'output',options.output);
out=dde_setupoutput('jac_res:reverse',options.output,out{:});
[J_rhs,res_rhs,ieq.rhs]=deal(out{:});
ieq.rhsvec=1:length(res_rhs);
neq=length(res_rhs);
%% combine res and J
res=[res_rhs;res_extra(:);zeros(length(options.step_cond),1)];
step_cond=repmat(point,numel(options.step_cond),1);
adj_cond=repmat(point,numel(options.step_cond),1);
pz=p_axpy(0,point,[]);
for i=numel(options.step_cond):-1:1
    step_cond(i)=dde_trim_point(feval(['dde_',point.kind,'_create'],options.step_cond(i)),...
        point);
    [dum,adj_cond(i)]=p_dot(pz,step_cond(i),'free_par_ind',free_par); %#ok<ASGLU>
end
J_cond=dde_x_from_point([cond_extra(:);adj_cond(:)],free_par);
J=[J_rhs;J_cond.'];
ieq.extra=neq+(1:length(res_extra));
neq=neq+length(res_extra);
ieq.step_cond=neq+(1:length(options.step_cond));
% neq=neq+(1:length(options.step_cond));
%% if extra parameters and extra column are provided, add them
if nextrapar>0 && extralen>0
    res(1:extralen)=res(1:extralen)+extracolumns*extrapar(:);
    J(1:extralen,end+(1:nextrapar))=extracolumns;
end
%% switch outputs if required
varargout=dde_setupoutput('jac_res',options.output,J,res,ieq);
end
