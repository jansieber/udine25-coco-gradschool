function dirderi=dde_dirderi_create(type,rhscell,varargin)
%% create cell of finite-difference derivative approximations 
% up to max_order using dde_dirderiv 
%
% argument type is 'tau' or 'rhs' to distinguish function formats
%
% argument rhscell contains cell with rhs function and already existing
% derivatives
%
% argument nf provides first dimension of output (only needed for rhs)
default={'maxorder',2,'hjac',@(ord)eps^(1/(2+ord)),'nf',[]};
options=dde_set_options(default,varargin,'pass_on');
derivbase=rhscell{end};
dirderi=rhscell(2:end);
if length(rhscell)==1 % not only rhs or sys_tau given (but lower derivative)
    derivbase=@(varargin)rhscell{1}(varargin{1:end-2});
    dirderi={};
end
ldirderi=length(dirderi);
nf=options.nf;
if isnumeric(nf)
    nf=@(x)nf;
end
for order=ldirderi+1:options.maxorder
    switch type
        case 'rhs'
            dirderi{order}=@(x,p,dx,dp)num_dirderiv(...
                @(xa,pa)derivbase(xa,pa,dx,dp),[2,2],order-ldirderi,...
                x,p,dx,dp,'nf',nf(x),'hjac',options.hjac(order));
        case 'tau'
            dirderi{order}=@(it,x,p,dx,dp)num_dirderiv(...
                @(xa,pa)derivbase(it,xa,pa,dx,dp),[2,2],order-ldirderi,...
                x,p,dx,dp,'nf',length(it),...
                'hjac',options.hjac(order));
    end
end
end
