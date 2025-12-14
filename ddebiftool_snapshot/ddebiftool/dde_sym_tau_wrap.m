function y=dde_sym_tau_wrap(fun,itau,order,xx,par,xdev,pdev)
%% Wrapper around automatically generated functions from symbolic differentiation
%
% $Id: dde_sym_tau_wrap.m 296 2018-09-24 21:36:56Z jansieber $
%% determine vectorized dimensions
% determine overall number of delays, which determines number of arguments
% for fun
xsize=size(xx);
nx=xsize(1);
%% find out which delay function to call and how to process its output
ndelays=fun('ntau')+1;
npar=fun('npar');
iscollected=fun('iscollected');
sys_tau_in=fun('sys_tau_in');
sys_tau_seq=fun('sys_tau_seq');
taunum=sys_tau_in(itau(1));
tau_exp=sys_tau_seq{taunum};
[dum,indout]=ismember(itau,tau_exp); %#ok<ASGLU>
nout=length(itau);
nexp=length(tau_exp);
%% prepare arguments
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
xapp=NaN(nx*(ndelays-min(itau)),nvec);
xx=cat(1,reshape(xx(:,1:min(itau),:),[nx*min(itau),nvec]),xapp);
directional_derivative=fun('directional_derivative');
if directional_derivative
    mfrep=1;
else
    maxorder=fun('maxorder');
    mfrep=numel(xdev)/numel(xx);
end
par=reshape(par,npar,[]);
xdev=cat(1,reshape(xdev(:,1:min(itau),:),[nx*min(itau)*mfrep,nvec]),xapp);
pdev=reshape(pdev,npar*mfrep,[]);
if ~directional_derivative && mfrep==1
    xdev=repmat(xdev,maxorder,1);
    pdev=repmat(pdev,[1,maxorder,1]);
    mfrep=maxorder;
end
%% convert arrays to cells/argument lists in first 2 or 3 dimensions
[args,dargs]=dde_cell_from_arrays(xx,xdev, par,pdev,...
    nx,ndelays,npar,mfrep,iscollected);
out=cell(1,nexp);
[out{:}]=fun('tau',taunum,order,nexp,args{:},dargs{:});
y=NaN(nout,nvec);
if ~iscollected
    out=cellfun(@(y){y.'},out);
end
for i=nout:-1:1
    y(i,:)=out{indout(i)};
end
y=reshape(y,[nout,vecdim]);
end
