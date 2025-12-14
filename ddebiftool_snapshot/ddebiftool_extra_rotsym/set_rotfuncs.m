%% set_rotfuncs - Fill in funcs structure for use with DDE-Biftool, rotational symmetric case
%%
function funcs=set_rotfuncs(varargin)
%% Named arguments
%
% * |'sys_rhs'| or |'wrap_rhs'| (mandatory), 
% * |'rotation'| (mandatory),
% * |'exp_rotation'| (mandatory),
% * |'expA_vectorized'| (opt,default false), if exp_rotation can be called
% with multiple phi.
% * |'sys_tau'| or |'wrap_tau'|, 
% * |'sys_dirdtau'| or |'wrap_dirdtau'|, 
% * |'sys_ntau'|, 
% * |'sys_cond'|,
%
% See rotsym_demo for example usage. Uses default finite differences for
% partial derivatives of order>1.
%
% $Id: set_rotfuncs.m 369 2019-08-27 00:07:02Z jansieber $
%

%% Process options
defaults={'rotation',[],'exp_rotation',[],'rot_tol',1e-8,'hjac',...
    @(ord)eps^(1/(2+ord)),...
    'expA_vectorized',false};
[options,pass_on]=dde_set_options(defaults,varargin,'pass_on');
rot_tol=options.rot_tol;
userfuncs=set_funcs(pass_on{:},'hjac',options.hjac);
%% check rotation matrices
if isempty(options.rotation)
    error('rotation matrix must be given as argument ''rotation''');
else
    userfuncs.rotation=options.rotation;
end
if isempty(options.exp_rotation)
    userfuncs.exp_rotation=@(phi)expm(userfuncs.rotation*phi);
    options.expA_vectorized=false;
else
    userfuncs.exp_rotation=options.exp_rotation;
end
userfuncs.vec_rotation=userfuncs.exp_rotation;
if ~options.expA_vectorized
    userfuncs.vec_rotation=...
        @(phi)vec_rotation(userfuncs.exp_rotation,phi);
end
% test compatibility
phi=linspace(0,2*pi,10);
expAphi=userfuncs.vec_rotation(phi);
for i=1:length(phi)
    err=expm(userfuncs.rotation*phi(i))-expAphi(:,:,i);
    if max(abs(err(:)))>rot_tol
        error('exp(rot*phi)~=exp_rotation(phi) for phi=%g',phi(i));
    end
end
err=userfuncs.exp_rotation(2*pi)-eye(size(userfuncs.rotation));
if max(abs(err(:)))>rot_tol
    error('exp(rot*2*pi)~=Id');
end
err=userfuncs.rotation+userfuncs.rotation.';
if max(abs(err(:)))>rot_tol
    error('rot~=-rot^T');
end
lhs_num=userfuncs.lhs_matrixfun(size(userfuncs.rotation,1));
if norm(lhs_num*userfuncs.rotation-userfuncs.rotation*lhs_num,'inf')>rot_tol
    error('lhs_matrix and rotation do not commute');
end
%% modify sys_rhs to include shift to rotating coordinates
rhs_args={userfuncs.rotation,userfuncs.vec_rotation,userfuncs,lhs_num};
rhs=@(xx,p)rot_rhs(rhs_args{:},0,xx,p);
%% set derivatives if provided by user
dirderi_args={};
if userfuncs.dirderi_provided()>0
    dirderi_args={'wrap_dirderi',...
        {@(xx,p,dx,dp)rot_rhs(rhs_args{:},1,xx,p,dx,dp)}};
end
funcs=set_funcs('wrap_rhs',rhs,dirderi_args{:},...
    'wrap_tau',userfuncs.wrap_tau,...
    'dtau_dir',userfuncs.dtau_dir,...
    'dtau_mf',userfuncs.dtau_mf,...
    'sys_ntau',dde_num_delays(userfuncs),...
    'lhs_matrix',lhs_num);
%% modify sys_cond to include fixing of rotational phase
rot_syscond=dde_sys_cond_create('name','rotationfix','fun',@(p,pref)rot_cond(p,userfuncs,pref),...
    'reference',true);
%% enable extracting underlying solution, removing omega
% and embedding (the reverse of extracting)
funcs.get_comp=@(pt,comp)extract_from_rot(pt,comp,size(userfuncs.rotation,1));
funcs.embed=@(p,component,template)rot_embed(p,component,template,size(userfuncs.rotation,1));
funcs.userfuncs=userfuncs;
funcs.vec_rotation=userfuncs.vec_rotation;
funcs.rotation=userfuncs.rotation;
embedded_conds=dde_embed_cond(funcs,funcs.userfuncs.sys_cond,'solution');
funcs.sys_cond=repmat(dde_sys_cond_create(),0,1);
funcs=dde_funcs_add_cond(funcs,embedded_conds,'name','rot_embedded');
funcs=dde_funcs_add_cond(funcs,rot_syscond);
end
%%
function rotarray=rot_embed(parray,component,template,dim)
if isfield(template,'x')
    sol='x';
else
    sol='profile';
end
        
rotarray=parray;
nx=size(template.(sol),1);
if strcmp(component,'omega')
    temp0=p_axpy(0,template,[]);
end
for i=numel(rotarray):-1:1
    rot=rotarray(i);
    switch component
        case {'solution','solution_for_stability'}
        rot.parameter(end+1)=0;
        x=rot.(sol);
        if nx>dim
            x(dim+1:nx,:)=0;
        end
        rot.(sol)=x;
    case 'omega'
        rot=temp0;
        rot.parameter(end)=p;
    end
    rotarray(i)=rot;
end
end
%%
function expAphi=vec_rotation(expA,phi)
for i=size(phi,2):-1:1
    expAphi(:,:,i)=expA(phi(i));
end
end
