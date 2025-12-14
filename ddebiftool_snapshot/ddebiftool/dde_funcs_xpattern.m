function xpattern=dde_funcs_xpattern(funcs,nx)
[ntau,max_rhs_xderiv]=dde_num_delays(funcs);
if nargin<2
    nx=size(funcs.lhs_matrixfun(),2);
end
xpattern=funcs.xpattern;
if ~isempty(xpattern)
    return
end
xpattern=dde_xpattern('incl_deriv',max_rhs_xderiv>0,'xdims',[nx,ntau+1,max_rhs_xderiv+1]);
end