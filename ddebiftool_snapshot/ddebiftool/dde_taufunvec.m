function [taus,ncalls]=dde_taufunvec(funcs,x,p,repeat,inc0)
%% compute all delays (incl 0) if x are already known
% used in dde_stst_delays and inside variational problems for psol
if nargin<5
    inc0=true;
end
if nargin<4
    repeat=false;
end
[xsize,psize]=deal(size(x),size(p));
vecdim=[xsize(3:end),1];
nvec=prod(vecdim);
x=reshape(x,xsize(1),xsize(2),nvec);
p=reshape(p,psize(1),psize(2),nvec);
xhistfun=@(it,t,deriv)x(:,it+1,:,deriv+1);
d_nr=1-double(inc0):dde_num_delays(funcs); 
[taus,xvals,ncalls]=...
    dde_fun_delays(funcs,xhistfun,size(x,3),p,'repeat',repeat,'d_nr',d_nr); %#ok<ASGLU>
taus=reshape(taus,[size(taus,1),vecdim]);
end
