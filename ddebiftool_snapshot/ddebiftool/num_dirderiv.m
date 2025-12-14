function [df,dflow]=num_dirderiv(fun,dims,order,varargin)
%% finite-difference directional derivatives approximations of finp of orbitrary order
%
% |num_dirderiv(f,na,j,x,...,dx,...)| returns approximately
%
% $$ \frac{\partial^j}{\partial h^j}f(x+h dx,...)\vert_{h=0}$
% 
%% Input
%
% * |finp| r.h.s, y=fun(x,...), which gets differentiated
% * |dims| number of dimensions of arguments of f
% (or|dims={outdim,argdims}|)
% * |order| order of requested derivative
% * |varargin{1:nargs}| arguments of f
% * |varargin{nargs+(1:nargs)}| deviation applied for derivative
%
%% Optional inputs
%
% * |hjac| (default eps^(1/(order+2))) initial scaling in spacing and output matrix
% * |fac| (1.99) ratio between successive deviations for high order
% derivatives
% * |nf| hint for output dimensions ignoring vectorization
%
%% Output
%
% * |df| approx derivative of order |order| at (xx,par) in direction
% (dx,dpar)
%
%% process options
if isnumeric(dims)
    argdims=dims;
    outdims=1;
else
    argdims=dims{2};
    outdims=dims{1};
end
argdims=argdims(:);
narg=length(argdims);
args=reshape(varargin(1:2*narg),narg,2);
x=reshape(varargin(1:narg),[],1);
dx=reshape(varargin(narg+(1:narg)),[],1);
fmt=cellfun(@(a,f){size(a,1:f)},x(:),num2cell(argdims));
sz1=size(x{1});
vecdim=[sz1(length(fmt{1})+1:end),1];
default={'hjac',eps^(1/(order+2)),'nf',[],'axpy',@(a,x,y)a.*x+y};
options=loc_set_opt(default,varargin(2*narg+1:end));
if all(cellfun(@isempty,args(:)))
    [df,dflow]=deal(zeros([options.nf,vecdim]));
    return
end
if ~isnumeric(options.hjac)
    options.hjac=options.hjac(order);
end
%% wrap for vectorization of fun
nvec=prod(vecdim);
xx=cellfun(@(a,f){reshape(a,[f,nvec])}, x,fmt);
dxx=cellfun(@(a,f){reshape(a,[f,nvec])},dx,fmt);
xdh=@(s)cellfun(@(x,d){options.axpy(s*options.hjac,d,x)},xx,dxx);
dflow=[];
switch order
    %% use simple formulas for low orders
    case 1
        [xp,xm]=deal(xdh(+1),xdh(-1));
        df=0.5*(fun(xp{:})-fun(xm{:}))/options.hjac;
    case 2
        [xp,xm]=deal(xdh(+1),xdh(-1));
        df=(fun(xp{:})+fun(xm{:})-2*fun(xx{:}))/options.hjac^2;
    otherwise
        nx=num2cell(cellfun(@prod,fmt));
        rhs=@(h)fdevh(fun,h,xx,dxx,fmt,nx,nvec,options.axpy);
        [df,dflow]=num_scalarderiv(rhs,order,...
            'h',options.hjac,varargin{2*narg+1:end});
        df=reshape(df,[],nvec);
        dflow=reshape(dflow,[],nvec);
end
if ~isempty(options.nf)
    nf=options.nf;
else 
    nf=size(df,1:outdims);
end
df=reshape(df,[nf,vecdim]);
if isempty(dflow)
    dflow=df;
end
end
function y=fdevh(f,h,xx,dx,fmt,nx,nvec,axpy)
nh=length(h);
xx=cellfun(@(a,f){repmat(reshape(a,f,nvec),[1,1,nh])},xx,nx);
dx=cellfun(@(a,f){repmat(reshape(a,f,nvec),[1,1,nh])},dx,nx);
hvec=cellfun(@(nxa){repmat(reshape(h,[1,1,nh]),[nxa,nvec,1])},nx);
xdev=cellfun(@(xa,ha,dxa){axpy(ha,dxa,xa)},xx,hvec,dx);
xdev=cellfun(@(a,f){reshape(a,[f,nvec*nh])},xdev,fmt);
y=f(xdev{:});
y=reshape(y,[size(y,1)*nvec,nh]);
end
%%
function opts=loc_set_opt(default,args)
opts=struct(default{:});
args=reshape(args,2,[]);
for i=1:size(args,2)
    if isfield(opts,args{1,i})
        opts.(args{1,i})=args{2,i};
    end
end
end
