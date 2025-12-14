function [taus,xtau_vals]=dde_coll_delays(funcs,pt,varargin)
%% compute all or requested delays in coll, including first delay=0
% at nt requested or mesh times t. taus is nt x d (d is number of delays
% incl 0). deriv_ord(k) is derivative of kth arg of sys_rhs, sys_tau etc,
% xtau_vals  is nx x d x nt x (maxderiv+1) array of all delayed values of
% profile and their derivatives
%%
ntaum1=dde_num_delays(funcs); 
pt_eval_choices=struct('psol',@psol_eval,'hcli',@hcli_eval);
default={'d_nr',0:ntaum1,'t',pt.mesh,'max_xderiv',1,'repeat',true,'pt_eval',pt.kind};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
tcheck=options.t(:)';
pt_eval=pt_eval_choices.(options.pt_eval);
xhistfun=@(it,negtau,deriv)pt_eval(pt,tcheck,it,negtau,0:options.max_xderiv,pass_on{:});
[taus,xtau_vals]=dde_fun_delays(funcs,xhistfun,...
    length(tcheck),pt.parameter,'repeat',options.repeat,'d_nr',options.d_nr,...
    'max_xderiv',options.max_xderiv);
end

%%
function y=psol_eval(pt,tbase,itau,negtau,orders,varargin)
ntauc=length(itau);
negtau=reshape(negtau,ntauc,[]);
t=tbase(ones(ntauc,1),:)+negtau/pt.period;
nderiv=length(orders);
xhfun=@(ord){dde_coll_wrap_eva(pt,t(:).',varargin{:},'diff',ord,'output','profile')};
y=cell2mat(arrayfun(xhfun,reshape(orders,1,1,[])));
y=reshape(y,size(y,1),ntauc,length(tbase),nderiv);
end
%%
function y=hcli_eval(pt,tbase,itau,negtau,orders,varargin)
[nx,ntauc,nb,nderiv]=deal(size(pt.profile,1),length(itau),size(pt.lambda_v,1),length(orders));
negtau=reshape(negtau,ntauc,[]);
t=tbase(ones(ntauc,1),:)+negtau/pt.period;
t=t(:).';
isneg=t<min(pt.mesh);
nneg=sum(isneg);
oneg=ones(1,nneg);
y=NaN(nx,ntauc*length(tbase),nderiv);
expDt=dde_expmt(pt.lambda_v,pt.period*t(isneg));
VexpDt=mshape(expDt,'rs',{nb,nb*nneg},'apply',@(M)pt.v*M,'rs',{nx,nb,nneg},...
    'p',[1,3,2],'rs',{nx*nneg,nb}); % (nx ntisneg) x nb
al=pt.alpha*pt.epsilon;
for i=1:length(orders)
    y(:,~isneg,i)=dde_coll_wrap_eva(pt,t(~isneg),varargin{:},...
        'diff',orders(i),'output','profile');
    Dpow=(pt.lambda_v*pt.period)^orders(i);
    y(:,isneg,i)=pt.x1(:,oneg).*double(orders(i)==0)+...
        reshape(VexpDt*Dpow*al,nx,nneg);
end
end