%% Codimension-2 normal form computation (wrapping for finite differences)
%% Input
%
% * funcs: problem functions
% * pt: point that is (approx) codim-2 point
% * type: type of normalform computation (used to check in pointtype_list
% for required order)
% * nmfm_compute: function of type newpoint=nmfm_type(funcs,point) where
% type is hoho, zeho, BT, cusp etc
%% Output
%
% * nf: point with normal form
% * nflow: if numerical finite differences are used then computation is
% done twice, once with higher order, once with lower order, this output is
% the result with lower order. Use the difference between nf and nflow to
% estimate the error
%
% $Id: nmfm_wrap_findiff.m 309 2018-10-28 19:02:42Z jansieber $
%%
function [nf,nflow]=nmfm_wrap_findiff(funcs,pt,type,nmfm_compute,varargin)
newpoint=nmfm_compute(funcs,pt,varargin{:});
nf=newpoint;
nflow=nf;
order=getfield(getfield(pointtype_list(),'order'),type);
if order>funcs.dirderi_provided()
    newpoint=nmfm_compute(funcs,pt,varargin{:},'output',2);
    nflow=newpoint;
end
end
