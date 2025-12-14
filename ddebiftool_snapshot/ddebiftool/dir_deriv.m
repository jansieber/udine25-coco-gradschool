function y=dir_deriv(fbase,fdirderivs,dims,order,varargin)
%% extend directional derivative to arbitrary order using finite differences
% fbase is function, fdirderivs
default={'hjac',@(ord)eps^(1/(ord+2)),'nf',[],'splitcomplex',false,'op',num_op()};
if isnumeric(dims)
    argdims=dims;
else
    argdims=dims{2};
end
argdims=argdims(:);
narg=length(argdims);
if nargin<=3 || ischar(order)
    [options,optcell]=loc_set_opt(default,[{order},varargin]); %#ok<ASGLU>
    y=@(ord,varargin)dir_deriv(fbase,fdirderivs,dims,ord,varargin{:},optcell{:});
    return
end
x=reshape(varargin(1:narg),[],1);
ndxend=narg+narg*double(order>0);
dx=reshape(varargin(narg+1:ndxend),[],1);
fmt=cellfun(@(a,f){size(a,1:f)},x(:),num2cell(argdims));
sz1=size(x{1});
vecdim=[sz1(length(fmt{1})+1:end),1];
[options,optcell]=loc_set_opt(default,varargin(ndxend+1:end));
% if all(cellfun(@isempty,[x(:);dx(:)]))
%     y=zeros([options.nf,vecdim]);
%     return
% end
if order==0
    y=fbase(x{:});
    return
end
if all(cellfun(@(x)all(x(:)==0),dx))
    y=zeros([options.nf,vecdim]);
    return
end
op=options.op;
hasimag=cellfun(@(a)any(~op.is0(op.imag(a(:)))),dx);
if order<=length(fdirderivs) && ~(options.splitcomplex&&any(hasimag))
    y=fdirderivs{order}(x{:},dx{:});
else
    %% if required split up complex dev arguments (dx) recursively for numerical approximation
    % to avoid complex base arguments
    if options.splitcomplex && any(hasimag)
        fun=@(varargin)dir_deriv(fbase,fdirderivs,dims,order,varargin{:},optcell{:});
        dxri=cat(2,cellfun(@(a){op.real(a)},dx(:)),...
            cellfun(@(a){op.imag(a)},dx(:)));
        y=mult_combine(fun,order,dims,op.matmul,x,dxri,'complex');
        return
    end
    excess_order=order-length(fdirderivs);    
    if ~isnumeric(options.hjac)
        options.hjac=options.hjac(excess_order);
    end
    if isempty(fdirderivs)
        y=num_dirderiv(fbase,dims,excess_order,...
            x{:},dx{:},'nf',options.nf,'hjac',options.hjac);
    else
        y=num_dirderiv(@(varargin)excess_dirderiv(fdirderivs{end},varargin(:),dx),...
            dims,excess_order,x{:},dx{:},'nf',options.nf,'hjac',options.hjac);
    end
end
end
%%
function y=excess_dirderiv(fdir,x,dx)
szx=cellfun(@(a){size(a)},x);
szdx=cellfun(@(a){size(a)},dx);
ratios=cellfun(@(xs,dxs){xs./(max(dxs,1))},szx,szdx);
dx=cellfun(@(a,r){repmat(a,r)},dx,ratios);
y=fdir(x{:},dx{:});
end
%%
function [opts,optcell]=loc_set_opt(default,args)
opts=struct(default{:});
args=reshape(args,2,[]);
for i=1:size(args,2)
    if isfield(opts,args{1,i})
        opts.(args{1,i})=args{2,i};
    end
end
optcell(1,:)=fieldnames(opts);
optcell(2,:)=cellfun(@(f){opts.(f)},optcell(1,:));
optcell=reshape(optcell,1,[]);
end
%%
%%
function op=num_op()
op.real=@(x)real(x);
op.imag=@(x)imag(x);
op.is0=@(x)x==0;
op.matmul=@(c,x)c*x;
end
