function [xpattern,xpatvec]=dde_xpattern(varargin)
default={'incl_deriv',false,'x',[],'xdims',[],'xpattern',[]};
options=dde_set_options(default,varargin,'pass_on');
[incl_deriv,x,xdims,xpattern]=deal(options.incl_deriv,options.x,options.xdims,options.xpattern);
assert(~isempty(x)||~isempty(xdims));
nxdims=2+double(incl_deriv);
if isempty(xdims)
    xdims=size(x,1:nxdims);
end
xlen=prod(xdims);
if isempty(xpattern)||any(isnan(xpattern(:)))
    xpatvec=1:xlen;
    xpatc=cell(nxdims,1);
    [xpatc{:}]=ind2sub(xdims,xpatvec);
    xpattern=cat(1,xpatc{:});
    return
end
if size(xpattern,1)==2 && nxdims==3
    xpattern=[xpattern;ones(1,size(xpattern,2))];
end
xpatc=num2cell(xpattern,2);
xpatvec=sub2ind(xdims,xpatc{:});
end