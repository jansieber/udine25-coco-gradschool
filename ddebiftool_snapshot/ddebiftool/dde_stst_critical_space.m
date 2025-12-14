function [out1,out2,out3]=dde_stst_critical_space(funcs,point,kind,omega,mth,varargin)
%% find critical subspace for eigenvalue 1i*omega (or 0)
default={'output','point','tolerance',1e-6,'v_scal',[],'extracolumns',[]};
[options,pass_on]=dde_set_options(default,varargin,'pass_on');
%% initialize point for extended system
pini=feval(['dde_',kind,'_create'],'point',point,'omega',omega);
pini.v=0*pini.x;
%% if defined find indices of free parameters that regularise problem
[free_par,ifree,iextfree,ip_extfree]=get_free_par(funcs);
%% define bifurcation problem and determine Jacobian Jext
if isfield(mth,'extracolumns')
    mth.extracolumns=[];
end
[f,x0]=p_correc_setup(funcs,pini,free_par,mth,'previous',pini,'output','J',pass_on{:});
[Jext,dum,ieq]=f(x0); %#ok<ASGLU>
%% find variational part of Jext, call it D
ind=dde_ind_from_point(pini,free_par);
[irow,n]=get_var_inds(ieq,{'rhs','dvar'});
icol=get_var_inds(ind,{'v'});
D=Jext(irow,[icol;reshape(ind.parameter(iextfree),[],1)]);
%% D should have 1d nullspace if number of constraints is correct
% assign this nullspace to variational state and parameters
[U,S,V]=svd(D); %#ok<ASGLU>
sv=diag(S);
v=V(:,end);
if strcmp(kind,'hopf')
    pini.v=v(1:n)+1i*v(n+1:2*n);
else
    pini.v=v(1:n);
end
pini.parameter(ip_extfree)=v(n+ifree);
%% Normalise variational part
if isempty(options.v_scal)
    resnorm=dde_stst_nullnorm_cond(pini,ip_extfree,mth,0);
else
    resnorm=options.v_scal(pini);
end
pini.v=pini.v/resnorm;
pini.parameter(ip_extfree)=pini.parameter(ip_extfree)/resnorm;
%% output, dependent on options
[out1,out2,out3]=deal(v,sv,pini);
switch options.output
    case {'sv','singular_values'}
        [out1,out2]=deal(sv,v);
    case {'dim','dimension'}
        sv=sum(sv<=options.tolerance);
        [out1,out2]=deal(sv,v);
    case {'pini'}
        out1=pini;
    case {'point'}
        [out1,out2,out3]=deal(pini,v,sv);
end
end
%%
function [free_par,ifree,iextfree,ip_extfree]=get_free_par(funcs)
[free_par,ifree,iextfree,ip_extfree]=deal(zeros(0,1));
if ~isfield(funcs,'ip') || ~isfield(funcs.ip,'nullparind') ||...
    isempty(funcs.ip.nullparind)
    return
end
free_par=funcs.ip.nullparind(:);
[dum,iextfree]=ismember(funcs.ip.nullparind(:,2),free_par); %#ok<ASGLU>
[dum,ifree]=ismember(funcs.ip.nullparind(:,1),free_par); %#ok<ASGLU>
ip_extfree=funcs.ip.nullparind(:,2);
end
%% extract rows for variational problem
function [irow,nvardim]=get_var_inds(ind,fields)
irow=ind;
for i=1:length(fields)
    irow=irow.(fields{i});
end
irow=irow(:);
nvardim=length(irow);
if isstruct(irow)
    fn=fieldnames(irow);
    nvardim=length(irow.(fn{1}));
    irow=cell2mat(cellfun(@(s){reshape(irow.(s),[],1)},fn(:)));
    irow=irow(:);
end
if isfield(ind,'extra')
    irow=[irow;ind.extra(:)];
end
end