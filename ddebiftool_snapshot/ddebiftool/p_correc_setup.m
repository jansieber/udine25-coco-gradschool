function [f,x0,data]=p_correc_setup(funcs,point,free_par,method,varargin)
%% set up f and initial guess x0,p reprocess and remesh point if required
%
% INPUT:
%
%  * funcs problem functions
%  * point initial point guess
%  * free_par free parameter numbers in N^d
%  * method: point method parameters
%  * step_cnd (default  empty) steplength condition(s) as point(s)
%  * remesh_flag (default 0) if zero or absent, do not adapt mesh; if one, always adapt
%  * previous (default empty) previously computed point along branch (for
%  example for periodic solutions or connecting orbits, to minimize phase
%  shift)
%
% OUTPUT:
%  * f function of type [r,J]=f(x) returning residual
%  * x0: vector, initial guess corresponding to point
%  * data: structure with remaining info (preprocessed free_par, method,
%  etc)
% 
% $Id: p_correc_setup.m 352 2019-06-19 16:03:03Z jansieber $
%%
default={'remesh_flag',false,'step_cnd',repmat(point,0,1),...
    'previous',repmat(point,0,1),'extracolumns',[],'row_minus_column',0};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% replace old stability and other info by empty fields if present:
point=dde_trim_point(feval(['dde_',point.kind,'_create'],point),point);
%% remesh if required
if isfield(method,'remesh') && ...
        method.remesh && ...
        options.remesh_flag>1 && ...
        mod(options.remesh_flag,method.adapt_mesh_before_correct)==0
    point=p_remesh(point);
    if ~isempty(options.previous)
        options.previous=dde_trim_point(p_remesh(options.previous,point),point);
    end
    for i=1:length(options.step_cnd)
        options.step_cnd(i)=dde_trim_point(p_remesh(options.step_cnd(i),point),point);
    end
end
%% turn step_cnd vectors into their adjoints
%for i=1:length(options.step_cnd)
%    [dum,stp_orth]=p_dot(p_axpy(0,options.step_cnd(i),[]),options.step_cnd(i)); %#ok<ASGLU>
%    options.step_cnd(i)=stp_orth;
%end
%% preprocess point, stepcond and functions
data=struct('funcs',funcs,'free_par',free_par,'method',method,...
    'step_cnd',options.step_cnd,'point',point,...
    'previous',options.previous,'extracolumns',options.extracolumns,'rJini',[]);
%% append extra columns if requested by method field
if isfield(method,'extracolumns') && ~isempty(method.extracolumns)
    [r,J0]=p_correc_rhs(data.funcs,data.method,data.point,data.free_par,...
        'pref',data.previous,'step_cond',data.step_cnd,...
        'extracolumns',[],pass_on{:},'output','res+J','x',[]);
    if isnumeric(method.extracolumns)&& method.extracolumns>0
        nulldim=method.extracolumns;
    elseif ischar(method.extracolumns) && strcmp(method.extracolumns,'auto')
        nulldim=max(size(J0,1)-size(J0,2)+options.row_minus_column,0);
    else
        nulldim=0;
    end
    if nulldim>0
        use_prev=~isempty(data.previous) && isfield(data.previous,'nvec') && ...
            isstruct(data.previous.nvec) && isfield(data.previous.nvec,'extracolumns') &&...
            size(data.previous.nvec.extracolumns,2)==nulldim && ...
            size(data.previous.nvec.extracolumns,1)>=size(J0,1);
        if use_prev
            extra_mat=cat(1,J0',data.previous.nvec.extracolumns(1:size(J0,1),:)');
            extra_rhs=cat(1,zeros(size(extra_mat,1)-nulldim,nulldim),eye(nulldim));
            if issparse(extra_mat)
                [L,U,P,Q]=lu(extra_mat);
                data.extracolumns=Q*(U\(L\(P*extra_rhs)));
            else
                data.extracolumns=extra_mat\extra_rhs;
            end
            [data.extracolumns,~]=qr(data.extracolumns,0);
        else
            data.extracolumns=dde_nullspaces_lr(J0','nulldim',nulldim);
        end
        J0=cat(2,J0,data.extracolumns);
        data.rJini={r,J0};
    end
end 
%% preprocess Newton iteration data if requested
if isfield(method,'preprocess') && ~isempty(method.preprocess)
    data=feval(method.preprocess,data);
end
%% prepare f and x0 for Newton-Raphson iterations:
f=@(x,varargin)p_correc_rhs(data.funcs,data.method,data.point,data.free_par,...
    'x',x,'pref',data.previous,'step_cond',data.step_cnd,...
    'extracolumns',data.extracolumns,pass_on{:},varargin{:});
if nargout>1
    x0=cat(1,dde_x_from_point(data.point,data.free_par),...
        zeros(size(data.extracolumns,2),1));
end
end
