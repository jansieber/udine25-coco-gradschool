function y=dde_sym_rhs_wrap(fun,order,xx,par,xdev,pdev)
%% Wrapper around automatically generated functions from symbolic differentiation
%
% $Id: dde_sym_rhs_wrap.m 315 2019-01-29 19:42:21Z jansieber $
%% determine vectorized dimensions
xsize=size(xx);
nx=xsize(1);
ndelays=fun('ntau')+1;
npar=fun('npar');
nf=fun('nf');
iscollected=fun('iscollected');
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
if isempty(xx)
    y=reshape([],[nx,vecdim]);
    return
end
xx=reshape(xx(:,1:ndelays,:),[nx*ndelays,nvec]);
par=reshape(par,npar,[]);
mfrep=numel(xdev)/numel(xx);
xdev=reshape(xdev(:,1:ndelays,:),[nx*ndelays*mfrep,nvec]);
pdev=reshape(pdev(1,1:npar*mfrep,:),npar*mfrep,[]);
%% if multi-derivative is called as directional derivative, replicate deviations
% if order>mfrep && ~fun('directional_derivative')
%     xdev=repmat(xdev,order/mfrep,1);
%     pdev=repmat(pdev,order/mfrep,1);
% end
%% convert arrays to cells/argument lists in first 2 or 3 dimensions
[args,dargs]=dde_cell_from_arrays(xx,xdev, par,pdev,...
    nx,ndelays,npar,mfrep,iscollected);
out=cell(1,nf);
if ~isempty(xdev)
    %% call derivative
    [out{:}]=fun('rhs',1,order,nf,args{:},dargs{:});
else
    %% call function
    [out{:}]=fun('rhs',1,order,nf,args{:});
end
%% The ith row gets either filled in or expanded
y=NaN(nf,nvec);
if ~iscollected
    out=cellfun(@(y){y.'},out);
end
for i=nf:-1:1
    y(i,:)=out{i};
end
y=reshape(y,[nf,vecdim]);
end

