function y=mult_from_dirderiv(fun,dims,varargin)
%% wrapper constructing derivatives from directional derivatives of arbitrary order
% input fun has form y=fun(k,dims,x,..,dx,...,opts) which permits
% na=length(dims) arguments x and for order>=1 na arguments dx for the
% directional deviation, returning directional derivative of order k.
% Argument ctrl has cell format: first element is number of output
% dimensions, 2nd element array of numbers of argument dimensions (needed
% for vectorization), 3rd (optional) element is structure op
% with fields 'plus','minus', matmul', which can be
% redefined if function argument are not numeric or symbolic: op.matmul
% multiplies ndev x ncases matrix of deviations with k x ndev matrix of
% real numbers.
%%
default={'op',num_op(),'groups',[]};
optstart=find([cellfun(@ischar,varargin),true],1,'first');
opts=loc_set_opt(default,reshape(varargin(optstart:end),2,[]));
if isnumeric(dims)
    dims={1,dims};
end
[outdims,argdims]=deal(dims{:});
argdims=argdims(:);
[args,fmt]=arg_array_expand(argdims,varargin{1:optstart-1});
baseargs=args(:,1);
dev=args(:,2:end);
vecdims=[size(baseargs{1},argdims(1)+1:ndims(baseargs{1})),1];
%% reshape arguments into fmt{i} x nvec
nvec=prod(vecdims);
order=size(dev,2);
if order==0
    baseargs=cellfun(@(a,f){reshape(a,[f,nvec,1])},baseargs,fmt(:,1));
    y=fun(baseargs{:});
else
    y=mult_combine(fun,order,dims,opts.op.matmul,baseargs,dev,opts.groups);
end
if length(vecdims)>1
    y=full(y);
end
y=reshape(y,[size(y,1:outdims),vecdims]);
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
%%
function op=num_op()
op.matmul=@(c,x)c*x;
end
