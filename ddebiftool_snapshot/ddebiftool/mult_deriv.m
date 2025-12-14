function y=mult_deriv(fun,varargin)
%% multidirectional and tensor derivatives based on directional derivatives
% as far as provided in fcell and finite differences above. Either call as
% mult_deriv(f,fdirderi,dims,...), where fdirderi is cell array of
% directional derivatives or mult_deriv(fdirderiv,...) where fdirderiv is
% function that can be called as fdirderiv(order,...) for arbitrary order
% (no dims argument in this call). The function fdirderiv will be
% constructed inside for the call with f,fdirderi,dims,... using dir_deriv.
% The alternative calling sequence is intended for cases where the user
% constructed a function fdirderiv using dir_deriv themselves.
%
% Other arguments:
%
% for mult_deriv(f,fdirderi,dims,x,...,dx1,...,dx2,...) or for
% mult_deriv(fdirderiv,dims,x,...,dx1,...,dx2,...) 
%
% dims is integer array specifying ndims of each argument x. All further
% dimensions of inputs x or dxj are interpreted as vectorised calls. Length
% of dims equals number of arguments for f, such that order of the
% derivative is determined by (length(varargin)-iarg)/length(dims)-1, where
% iarg=3 for call type 1 and irg=2 for call type 2.
%
% If called without x and dx arguments a function handle for df is
% returned that can be used as df(x{:},dx{:}).
%%
if iscell(varargin{1})&&(isempty(varargin{1})||~isnumeric(varargin{1}{1})) % calling sequence mult_deriv(fun,fdirderi,dims,...)
    dfdir_inp=false;
    iarg=3;
    inp={fun,varargin{1}};
    dims=varargin{2};
else % calling sequence mult_deriv(fdirderiv,dims,....)
    dfdir_inp=true;
    iarg=2;
    dims=varargin{1};
end
%% optional inputs are hjac and nf
if isnumeric(dims)
    argdims=dims;
else
    argdims=dims{2};
end
argdims=argdims(:);
optstart=find([cellfun(@ischar,varargin),true],1,'first');
narg=length(argdims);
opts=reshape(varargin(optstart:end),2,[]);
%% in both cases order is first arg of dfdir
if ~dfdir_inp
    dfdir=dir_deriv(inp{:},dims,opts{:});
else
    dfdir=fun;
end
if optstart-1>=iarg
    y=loc_mult_deriv(dfdir,dims,narg,opts,varargin{iarg:optstart-1});
else
    y=@(varargin)loc_mult_deriv(dfdir,dims,narg,opts,varargin{:});
end
end
%%
function y=loc_mult_deriv(dfdir,dims,narg,opts,varargin)
args=reshape(varargin,narg,[]);
[base,dev]=deal(args(:,1),args(:,2:end));
order=size(dev,2);
y=mult_from_dirderiv(@(varargin)dfdir(order,varargin{:}),dims,base{:},dev{:},opts{:});
end