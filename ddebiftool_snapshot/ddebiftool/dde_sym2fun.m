function out=dde_sym2fun(fcnrhs,fun)
%% extract rhs  or delays or derivatives from symfuncs file
% fcnrhs is filename or function hadnle genreated by dde_sym2funcs fun
% (default 'rhs') indicates if delay or rhs should be returned
if nargin<2
    fun='rhs';
end
if ischar(fcnrhs)
    f=str2func(fcnrhs);
else
    f=fcnrhs;
end
dir_deriv=f('directional_derivative');
npar=f('npar');
switch fun
    case {'rhs','sys_rhs'}
        if dir_deriv
            out=@(x,p)dde_sym_rhs_wrap(f,0,x,p(1,1:npar,:),x,p(1,1:npar,:));
        else
            out=@(x,p)dde_sym_rhs_wrap(f,0,x,p(1,1:npar,:),{},{});
        end
    case {'dirderi','deri','rhs_dirderi'}
        %assert(dir_deriv);
        maxorder=f('maxorder');
        out=arrayfun(@(order){...
            @(x,p,dx,dp)dde_sym_rhs_wrap(f,order,x,p(1,1:npar,:),dx,dp(1,1:npar,:))},...
            1:maxorder);
    case 'mfderi'
        assert(~dir_deriv);
        out=@(xx,p,varargin)wrap_mfderi(f,xx,p,varargin);
    case {'delay','tau','sys_tau'}
        assert(dir_deriv);
        out=@(itau,x,p)dde_sym_tau_wrap(f,itau,0,x,p(1,1:npar,:),x,p(1,1:npar,:));
    case {'dirdtau','tauderi','tau_dirderi'}
        assert(dir_deriv);
        maxorder=f('maxorder');
        out=arrayfun(@(order){...
            @(itau,x,p,dx,dp)dde_sym_tau_wrap(...
            f,itau,order,x,p(1,1:npar,:),dx,dp(1,1:npar,:))},...
            1:maxorder);
    case {'ntau','sys_ntau'}
        assert(dir_deriv);
        out=@()ntau;
    otherwise
        error('dde_sym2fun:undefined',...
            'dde_sym2fun: function %s undefined',fun)
end
end