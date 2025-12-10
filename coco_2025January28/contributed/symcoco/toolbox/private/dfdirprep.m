%% wrapper for directional derivatives, expanding directions not provided
function y=dfdirprep(fun,fmt,varargin)
fargs=varargin(1:fmt.nargs);
devargs=varargin(fmt.nargs+1:end);
if isempty(devargs)
    devargs=repmat({},1,fmt.nargs);
end
[u0,udev,dims,argdim]=arg_expand_directional(fmt,fargs{1:fmt.nargs},devargs{:});
y=dfdir_wrap(fun,fmt,u0,udev,dims,argdim);
end
