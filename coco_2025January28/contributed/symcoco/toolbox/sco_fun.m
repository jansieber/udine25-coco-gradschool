function fout=sco_fun(fun_inp,argnames,varargin)
%% Return function or its derivative when directional derivatives or symbolic formulas are provided by user
%% 
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
%%
hdev0=eps.^(1./[3,4,6]);
default={'vector',true(1,length(argnames)),'hdev',hdev0,'debug',false};
options=sco_set_options(default,varargin,'pass_on');
fun=fun_inp;
if ~iscell(fun_inp)
    fun={fun_inp};
end
if length(options.hdev)<3
    options.hdev(end+1:3)=hdev0(length(options.hdev)+1:3);
end
if isscalar(fun)
    for i=1:3
        fun{i+1}=@(varargin)dirderiv(i,false,fun{1},options.hdev(i),...
        varargin(1:length(varargin)/2),varargin(length(varargin)/2+1:end));
    end
else
    nf=length(fun);
    for i=1:3
        fun{i+nf}=@(varargin)dirderiv(i,true,fun{nf},options.hdev(i),...
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
function y=dirderiv(ord,ext,f,h,args,devs)
cf=zeros(3,7);
mid=4;
cf(1,mid+([-1,1]))=[-1,1]/2;
cf(2,mid+([-1,0,1]))=[1,-2,1];
cf(3,mid+([-3:-1,1:3]))=[1,-8,13,-13,8,-1]/8;
evh=-3:3;
evh=evh(cf(ord,:)~=0);
cf_used=cf(ord,cf(ord,:)~=0);
y=0;
if ext
    df=devs;
else
    df={};
end
for i=1:length(evh)
    a=cellfun(@(x,d){x+evh(i)*h*d},args,devs);    
    y=y+cf_used(i)*f(a{:},df{:});
end
y=y/h^ord;
end
