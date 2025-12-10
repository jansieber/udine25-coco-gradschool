%% wrapper for named partial derivatives, converting them to directional derivatives
function y=dfnamed(fun,args,ndirs,fmt,varargin)
nargs=size(args,1);
args(ndirs)=varargin(nargs+(1:length(ndirs)));
devargs=arrayfun(@(i)cat(2,args(i,:)),1:nargs,'UniformOutput',false);
[u0,udev,dims,argdim]=arg_expand_mixed_derivative(fmt,{varargin{1:nargs},devargs{:}}); %#ok<CCAT>
y=dfdir_wrap(fun,fmt,u0,udev,dims,argdim);
end
