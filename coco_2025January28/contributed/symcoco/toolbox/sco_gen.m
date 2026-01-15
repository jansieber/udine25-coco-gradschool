function fout=sco_gen(fun,varargin)
%% Return function or its derivative stored in symbolical toolbox created fun, 
% or its derivative when directional derivatives are provided by user
%
% test is performed by call of fun('maxorder'): if it fails, we have
% user-provided function
%
%% For user-provided and finite-difference usage
% *Inputs:*
%
% * |fun|: function name or cell array
% * |argnames|: names for arguments, indicating number of arguments.
%
% *Optional/named Inputs:*
%
% * |debug| (default: false): assertions (which cost some time) are switched on
% * |vector| (default:true(1,length(argnames)): treat input as vector (when
% creating derivative tensors)
% * |hdev| (default (eps.^(1/[3,4,6]~[6e-6,1e-4,2e-3] deviations to be used
% ("h") for finite differnces if derivatives are not provided
%
% *Outputs:* function handle |fout|,which can be called with the number and
% format of arguments indicated by |args|. At the moment only functions
% with column vector inputs and a single column vector output are
% supported. All functions are vectorized, such that, e.g.,
%
% * After |F=sco_fun(fun,{'x','p'})|, |f=F('')|, |y=f(x,p)| has output |y| with
% |size(x,2)==size(p,2)| columns. Single-column expansion is enabled.
% 
% * After |df=F('x')|, |dy=df(x,p)| has output |dy| with
% |size(dy,2)==size(x,1)|, |size(dy,3)==size(x,2)==size(p,2)|.
%
% After |df2=F({'x','p'}|, |dy=df2(x,p)| has output |dy| with
% |size(dy,2)==size(x,1)|, |size(dy,3)==size(p,1)|,
% |size(dy,4)==size(x,2)==size(p,2)|.
%
% If ones of the input arguments has single column it will get expanded by repmat.
%
% Directional derivatives: if a '*' is present in any of the arguments of
% sco_gen then a directional derivative in this direction is returned. For
% example, after
% dfxvp=F({'x*v','p'}),dy=dfxvp(x,p,v) output dy has
% |size(dy,2)==size(p,1), |size(y,3)==max(size(x,2),size(p,2))|. It equals
% $\partial_{xp}f(x,p)v(.)$ where v is the deviation wrt x.
%
%% For symbolic toolbox usage
% *Inputs:*
%
% * |fun|: function name (or filename) of function created by sco_sym2funcs
% * |name|: if not present, |f| is returned, if char array or single cell with
% character, first derivative of |f| with respect to this argument is
% returned, if name is cell of length two, second derivative of |f| is
% returned. If name is numeric integer k then the directional derivative of
% order k is returned.
% * |debug|: if present assertions (which cost some time) are switched on
%
% If argument |name| is not present,
% |F=@(varargin)sco_gen(fun,varargin{:})| is returned, which can be used as
% an abbreviated call to |sco_gen|. E.g., |F('x')| is the same as
% |sco_gen(fun,'x')| after this initial call.
%
% number, format and names of arguments of functions can be checked with
% call |args=fun('argrange');|, which returns a struct |args|.
%
% *Outputs:* function handle |fout|,which can be called with the number and
% format of arguments indicated by |args|. At the moment only functions
% with column vector inputs and a single column vector output are
% supported. All functions are vectorized, such that, e.g.,
%
% * After |f=sco_gen(fun,'')|, |y=f(x,p)| has output |y| with
% |size(x,2)==size(p,2)| columns. Single-column expansion is enabled.
% 
% * After |df=sco_gen(fun,'x')|, |dy=df(x,p)| has output |dy| with
% |size(dy,2)==size(x,1)|, |size(dy,3)==size(x,2)==size(p,2)|.
%
% After |df2=sco_gen(fun,{'x','p'}|, |dy=df2(x,p)| has output |dy| with
% |size(dy,2)==size(x,1)|, |size(dy,3)==size(p,1)|,
% |size(dy,4)==size(x,2)==size(p,2)|.
%
% If ones of the input arguments has single column it will get expanded by repmat.
%
% Directional derivatives: if a '*' is present in any of the arguments of
% sco_gen then a directional derivative in this direction is returned. For
% example, after
% dfxvp=sco_gen(fun,{'x*v','p'}),dy=dfxvp(x,p,v) output dy has
% |size(dy,2)==size(p,1), |size(y,3)==max(size(x,2),size(p,2))|. It equals
% $\partial_{xp}f(x,p)v(.)$ where v is the deviation wrt x.
%%
try
    maxorder=fun('maxorder'); %#ok<NASGU>
    fout=sco_sym(fun,varargin{:});
catch
    fout=sco_fun(fun,varargin{:});
end
end
