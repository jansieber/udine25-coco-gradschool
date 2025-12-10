function fout=sco_fun(fun_inp,argnames,varargin)
%% Return function or its derivative when directional derivatives or symbolic formulas are provided by user
%% 
% *Inputs:*
%
% * |fun|: function name or cell array
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
default={'vector',true(1,length(argnames)),'output','fun',...
    'hdev',eps.^(1./[3,4]),'debug',false};
options=sco_set_options(default,varargin,'pass_on');
fun=fun_inp;
if ~iscell(fun_inp)
    fun={fun_inp};
end
if isscalar(fun)
    fun{2}=@(varargin)dirderiv1(fun{1},options.hdev(1),...
        varargin(1:length(varargin)/2),varargin(length(varargin)/2+1:end));
    fun{3}=@(varargin)dirderiv2(fun{1},options.hdev(2),...
        varargin(1:length(varargin)/2),varargin(length(varargin)/2+1:end));
else
    nf=length(fun);
    for i=nf+(1:2)
        fun{i}=@(varargin)dirderiv1ext(fun{end},options.hdev(1),...
            varargin(1:length(varargin)/2),varargin(length(varargin)/2+1:end));
    end
end    
fout=@(name)f_generate(fun,options,argnames,name);
end
%%
function fout=f_generate(fun,options,argnames,name)
[maxorder,nargs]=deal(length(fun)-1,length(argnames));
fmt_fun=@(order)arg_fmt(order,options.debug,nargs,options);
if (ischar(name) && strcmp(name,'')) || ... % F('') is function
        (iscell(name) && isempty(name)) || ...  % or F({})
        (iscell(name) && ischar(name{1}) && isempty(name{1})) % or F({''})
    fmt=fmt_fun(0);
    fout=@(varargin)dfdirprep(fun,fmt,varargin{:});
elseif isnumeric(name) && ismember(name,0:maxorder) % F(1), F(2),...
    fmt=fmt_fun(name);
    fout=@(varargin)dfdirprep(fun,fmt,varargin{:});
elseif ischar(name) || (iscell(name)&&ischar(name{1})) % F('x'), F({'x'}), F({'w','p'}) 
    [args,ndirs]=assemble_directions(argnames,name);
    fmt=fmt_fun(size(args,2));
    fout=@(varargin)dfnamed(fun,args,ndirs,fmt,varargin{:});
elseif iscell(name) && isscalar(name) && isnumeric(name{1})
    fmt=arg_fmt(fun,name{1},debug);
    fout=@(varargin)dfdir_wI(fun,fmt,varargin{:});
else    
    error('sco_fun:arg',['sco_fun: second argument ''name'' not',...
        'recognized, only derivatives up to order %d implemented'],maxorder);
end
end
%% Wrapper around automatically generated functions from symbolic differentiation
% converts numerical arrays into lists of scalar/vectorized arguments, as
% this is what the output from the symbolic toolbox produces.
function y=fuwrap(fun,order,argdim,u,du)
%% determine vectorized dimensions
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
%% 
function fmt=arg_fmt(deg,debug,nargs,options)
fmt.deg=deg;
fmt.nargs=nargs;
fmt.isvec=options.vector;
fmt.debug=debug;
fmt.fuwrap=@fuwrap;
end
%%
function y=dirderiv1(f,h,args,devs)
ap=cellfun(@(x,d){x+h*d},args,devs);
am=cellfun(@(x,d){x-h*d},args,devs);
y=(f(ap{:})-f(am{:}))/(2*h);
end
%%
function y=dirderiv1ext(f,h,args,devs)
ap=cellfun(@(x,d){x+h*d},args,devs);
am=cellfun(@(x,d){x-h*d},args,devs);
y=(f(ap{:},devs{:})-f(am{:},devs{:}))/(2*h);
end
%%
function y=dirderiv2(f,h,args,devs)
ap=cellfun(@(x,d){x+h*d},args,devs);
am=cellfun(@(x,d){x-h*d},args,devs);
y=(f(ap{:})+f(am{:})-2*f(args{:}))/h^2;
end
%%
function y=dirderiv2ext(f,h,args,devs)
ap=cellfun(@(x,d){x+h*d},args,devs);
am=cellfun(@(x,d){x-h*d},args,devs);
y=(f(ap{:},devs{:})+f(am{:},devs{:})-2*f(args{:},devs{:}))/h^2;
end
