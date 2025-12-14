function [taus,xtau_vals,ncalls]=...
    dde_fun_delays(funcs,xhistfun,argdim,par,varargin)
%% obtain values at -tau(i) for hist(t)
% argdim is number (w shape) of points where delay is to be evaluated. So,
% xhist will be called as xhist(itau,t,deriv_orders) where t is of size 1 x
% nt, itau indicates which delay number is being requested (starting count
% from 0), and nt=prod(argdim). It should return nx x nt x
% numel(deriv_orders) values. par has shape 1 x npar x {argdim or 1}
ntaum1=dde_num_delays(funcs); 
default={'d_nr',0:ntaum1,'max_xderiv',0};
options=dde_set_options(default,varargin,'pass_on');
ndel_out=length(options.d_nr);
max_xderiv=options.max_xderiv;
cellargs=arg_array_expand([2,1],par,NaN([1,argdim]));
par_vec=cellargs{1};
nvec=size(par_vec,3);
if ~funcs.tp_del % constant delay case
    tau=par_vec(1,funcs.sys_tau(),:);             % delay values
    tau=reshape(tau,size(tau,2),size(tau,3));
    taus=[zeros(1,size(tau,2));tau];   % delays (incl tau=0) x evaluations
    taus=taus(options.d_nr+1,:);
    if nargout>=2
        xh=xhistfun(options.d_nr,-reshape(taus,1,[]),0:max_xderiv);
        xtau_vals=reshape(xh,[size(xh,1),ndel_out,argdim,max_xderiv+1]);
        xtau_vals=xtau_reshape(xtau_vals,max_xderiv);
    end
    taus=reshape(taus,[ndel_out,argdim]);
    ncalls=1;
    return
end
taucalls=[{1},cellfun(@(i){i+1},funcs.sys_tau_seq)];
itau_max=max(options.d_nr+1);
taus=[zeros(1,nvec);NaN(ntaum1,nvec)];
ncalls=find(cellfun(@(itau)ismember(itau_max,itau),taucalls),1,'first')-1;
ndel_out=taucalls{ncalls+1}(end);
xh=xhistfun(0,-taus(1,:),0:max_xderiv);
nx=size(xh,1);
xtau_vals=cat(2,reshape(xh,nx,1,nvec,max_xderiv+1),...
    NaN(nx,ndel_out-1,nvec,max_xderiv+1));
for nc_i=2:ncalls+1
    taurg=1:max(taucalls{nc_i-1});
    yeval=xtau_vals(:,taurg,:,1);
    taus(taucalls{nc_i},:)=funcs.wrap_tau(taucalls{nc_i}-1,yeval,par_vec);
    t_i=taucalls{nc_i};
    ntau_i=length(t_i);
    xh=xhistfun(t_i-1,-reshape(taus(t_i,:),1,[]),0:max_xderiv);
    xtau_vals(:,t_i,:,:)=reshape(xh,[nx,ntau_i,nvec,max_xderiv+1]);
end
taus=taus(options.d_nr+1,:);% select requested delays
xtau_vals=xtau_reshape(xtau_vals(:,options.d_nr+1,:,:),max_xderiv);
end
%%
function xtaur=xtau_reshape(xtau,max_xderiv)
xtaur=xtau;
if max_xderiv==0
    return
end
xtaur=permute(xtaur,[1,2,4,3]);
end
