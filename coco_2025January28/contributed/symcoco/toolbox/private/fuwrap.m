%% Wrapper around automatically generated functions from symbolic differentiation
% converts numerical arrays into lists of scalar/vectorized arguments, as
% this is what the output from the symbolic toolbox produces.
function y=fuwrap(fun,order,argdim,u,du)

uargs=mat2cell(u,argdim,size(u,2));
if order==0
    y=fun{1}(uargs{:});
    return
end
if nargin<=4
    du=zeros(size(u));
end
duargs=mat2cell(du,argdim,size(u,2));
y=fun{order+1}(uargs{:},duargs{:});
end
