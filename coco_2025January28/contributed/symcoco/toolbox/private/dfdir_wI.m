%% wrapper for directional derivatives, expanding directions indicated by 'I'
function y=dfdir_wI(fun,fmt,varargin)
devargs=varargin(fmt.nargs+1:end);
[u0,udev,dims,argdim]=arg_expand_mixed_derivative(fmt,{varargin{1:fmt.nargs},devargs{:}}); %#ok<CCAT>
y=dfdir_wrap(fun,fmt,u0,udev,dims,argdim);
end
