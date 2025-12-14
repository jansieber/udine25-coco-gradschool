function J=dde_gen_deriv(dffunc,xx,par,nx,np,v)
%% compose derivatives of r.h.s wrt state and parameter from directional derivatives
%
% function J=dde_gen_deriv(dirderi,xx,par,nx,np,v)
%
%% INPUT:
%   dffunc: drhs_mf constructed from sys_dirderi:
%   dffunc(x,p,dx1,dp1,dx2,dp2) permitted, with expansion
%	xx state variable and delayed state variables columnwise
%	par list of parameter values
%	nx empty or list of requested state-derivatives (numbers of delay or zero)
%	np empty or list of requested parameter-derivatives
%	v matrix to multiply result with for 2nd xx derivative
%% OUTPUT:
%	J result of derivatives on righthandside (multiplied with v for 2nd
%	derivative of x) in old sys_deri format
%
%
%%
xsize=size(xx);
n=xsize(1);
ndelays=xsize(2);
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
xx=reshape(xx,[n,ndelays,1,nvec]);
par=reshape(par,1,size(par,2),1,[]);
if length(nx)==1 && isempty(np) && isempty(v)
    %% first order derivatives of the state:
    dx=dx_eye(n,ndelays,nx);
    df=dffunc(xx,par,dx,0);
    J=reshape(df,[size(df,1),n,vecdim]);
elseif isempty(nx) && length(np)==1 && isempty(v)
    %% first order derivative with respect to parameter
    J=dffunc(xx,par,0,{1,{np}});
    J=reshape(J,[size(J,1),1,vecdim]);
elseif length(nx)==2 && isempty(np) && ~isempty(v)
    %% second order state derivatives, applied to deviation
    dx=dx_eye(n,ndelays,nx(1));
    vext=zeros(n,n,1,nvec);
    v=reshape(v,n,1,1,[]);
    vext(:,nx(2)+1,1,:)=repmat(v,1,1,1,nvec/size(v,3));
    df=dffunc(xx,par,dx,0,vext,0);
    J=reshape(df,size(df,1),n,vecdim);
elseif length(nx)==1 && length(np)==1 && isempty(v)
    %% mixed state-parameter derivatives:
    xx=reshape(xx,[n,ndelays,1,1,nvec]);
    par=reshape(par,1,size(par,2),1,1,[]);
    dx=dx_eye(n,ndelays,nx);
    df=dffunc(xx,par,dx,0,0,{2,{np}});
    J=reshape(df,[size(df,1),n,vecdim]);
elseif isempty(nx) && length(np)==2
    %% 2nd-order derivative with respect to parameter
    xx=reshape(xx,[n,ndelays,1,1,nvec]);
    par=reshape(par,1,size(par,2),1,1,[]);
    df=dffunc(xx,par,0,{1,{np(1)}},0,{2,{np(2)}});
    J=reshape(df,[size(df,1),1,vecdim]);
else
    %% not implemented
    sv=num2str(size(v));
    error('dde_gen_deriv:notimplented',...
        ['dde_gen_deriv: requested derivative nx=%d, np=%d, size(v)=(%s)',...
        ' does not exist!'],nx,np,sv);
end
end
%%
function dx=dx_eye(n,ndelays,nxdev)
dx=zeros(n,n,ndelays);
dx(:,nxdev+1,:)=eye(n);
end