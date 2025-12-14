function varargout=arg_flatten(fmt,varargin)
%% convert arg list with uniform shapes of vectorization to argshape x nvec
% returns vecdim (shapes of vectrized dimensions) for later conversion back
if isempty(varargin)
    varargout={[],0};
    return
end
if iscell(varargin{1})
    args=varargin{1};
    cellout=true;
else
    args=varargin;
    cellout=false;
end
vecdim=[size(args{1},fmt(1)+1:ndims(args{1})),1];
nvec=prod(vecdim);
argsz=size(args);
args=reshape(args,length(fmt),[]);
fmt=repmat(fmt(:).',1,size(args,2));
args=arrayfun(@(i){reshape(args{i},[size(args{i},1:fmt(i)),nvec])},1:numel(args));
args=reshape(args,argsz);
if cellout
    varargout={args,vecdim,nvec};
else
    varargout=[args(:)',{vecdim},{nvec}];
end
end