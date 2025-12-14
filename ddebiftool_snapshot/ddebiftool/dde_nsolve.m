function [x,success]=dde_nsolve(f,x0,method,varargin)
%% Basic Newton iteration used in p_correc, solving f(x)=0 from x0, 
% independent of DDE-Biftool details such that it can be replaced by other
% routines or re-used outside of p_correc.
%
%  (c) DDE-BIFTOOL v. 3.2a(309), 28/10/2018
%%
default={'rJini',[]}; % sometimes the setup computes already the first s and J
options=dde_set_options(default,varargin,'pass_on');
max_iter=method.newton_max_iterations; % max number of newton iterations
nmon_iter=method.newton_nmon_iterations; % max number of nonmonotone iterations
conv_r=method.halting_accuracy; % required accuracy for stopping
min_r=method.minimal_accuracy; % minimal accuracy requred for convergence
print_r=method.print_residual_info; % print residual evolution
if isfield(method,'skip_jacobian')
    skip_jacobian=method.skip_jacobian;
else
    skip_jacobian=0;
end
x=x0;
[dx,r1,J]=deal(NaN(size(x)),NaN,NaN);
for i=1:max_iter
    try
        [res,J]=eval_fun(f,x,J,i,r1,options.rJini,skip_jacobian);
    catch ME
        success=false;
        if print_r
            fprintf(',\nillegal argument of nonlinearity\n%s\n',ME.message);
        end
        return
    end
    if ~isfield(method,'jacobian_nonsquare') || ~method.jacobian_nonsquare
        if size(J,1)~=size(J,2)
            warning('p_correc:nonsquare','dde_nsolve warning: use of nonsquare Jacobian.');
        end
    end
    dx_old=dx;
    dx=linsolve(J,-res);
    norm_res=norm(res,'inf');
    norm_dx=norm(dx,'inf');
    if any(isnan(dx))&& norm_res<min_r
        dx=dx_old;
        break;
    end
    x=x+dx;
    if print_r
        if i>1
            fprintf('\n');
        end
        fprintf('it=%3d, res=%8.4e, dx=%8.4e',i, norm_res, norm_dx);
    end
    % check for convergence
    if i==1
        r1=norm_res;
        r0=r1/2;
    else
        r0=r1;
        r1=norm_res;
    end
    if r1<=conv_r || (i-1>nmon_iter && r1>=r0) || any(isnan(r1)) || any(isnan(dx))
        break
    end 
end
% compute final residual
n_res=norm(res,'inf');
success=(n_res<=min_r) && all(~isnan(dx));
if print_r
    fprintf(', success=%2d\n',success);
end
if isfield(method,'print_jacobian_condition') && ...
        method.print_jacobian_condition && exist('J','var')
    fprintf('norm(J^(-1))=%g\n',1/min(svd(J)));
end
end
%%
function x=linsolve(A,y)
if issparse(A)
    [L,U,P,Q,R]=lu(A);
    x=Q*(U\(L\(P*(R\y))));
else
    x=A\y;
end
end
%%
function [res,J]=eval_fun(f,x,J,i,r1,rJini,skip_jacobian)
if i==1
    if ~isempty(rJini)
        [res,J]=deal(rJini{:});
    else
        [res,J]=f(x);
    end
else
    if skip_jacobian>0
        res=f(x,'output','res');
        norm_res=norm(res,'inf');
        if norm_res>=r1*skip_jacobian
            [res,J]=f(x);
        end
    else
        [res,J]=f(x);
    end
end
end