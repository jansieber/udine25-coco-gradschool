function [g,ng,name]=dde_sym_dist_delay(ginp,varargin)
%% wrapper: preprocess if g provided by dde_symfuncs
% either provide integrand or cell array of integrand and directional
% derivatives, or handle/function file name of dde_sym2funcs filename if
% integrand is provided, option 'nf' (dimension of output) is mandatory
default={{'gt','g_of_t_x','t_is_arg'},true,'name','',{'nf','ng'},NaN};
options=dde_set_options(default,varargin,'pass_on');
g_is_sym=false;
if ~iscell(ginp)
    g_is_sym=true;
    try
        gfun=dde_sym2fun(ginp,'rhs');
    catch
        g_is_sym=false;
        ginp={ginp};
    end
end
if g_is_sym
    gderi=dde_sym2fun(ginp,'dirderi');
    if isempty(options.name)
        if ischar(ginp)
            name=ginp;
        else
            name=func2str(ginp);
        end
    end
    ng=ginp('nf');
else
    gfun=ginp{1};
    gderi=ginp(2:end);
    if isempty(options.name)
        name=func2str(gfun);
        if ~isvarname(name)
            name='dist_delay';
        end
    end
    ng=options.nf;
end    
%% if the kernel depends on s, it should be the first component of the x 
% argument in the symbolic function
g=cell(1,1+length(gderi));
if options.gt
    if g_is_sym
        g{1}=@(s,x,p)gfun(cat(1,s,x),p);
        for i=1:length(gderi)
            g{i+1}=@(s,x,p,ds,dx,dp)gderi{i}(cat(1,s,x),p,cat(1,ds,dx),dp);
        end
    else
        g=[{gfun},gderi];
    end
else
    g{1}=@(s,x,p)gfun(x,p);
    for i=1:length(gderi)
        g{i+1}=@(s,x,p,ds,dx,dp)gderi{i}(x,p,dx,dp);
    end
end
end