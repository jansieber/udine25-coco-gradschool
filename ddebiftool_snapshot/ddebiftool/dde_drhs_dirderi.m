%% create full jacobians or rhs from directional derivatives
function drhs=dde_drhs_dirderi(funcs,x,par,varargin)
default={'output','x','derivs',4,'free_par',[]};
options=dde_set_options(default,varargin,'pass_on');
[ntau,deriv_order]=dde_num_delays(funcs);
vecdim=[size(x,setdiff(3:ndims(x),options.derivs)),1];
nvec=prod(vecdim);
[npar,npv]=size(par,[2,3]);
p=par(1,:,ones(1,nvec/npv));
ntaux=(ntau+1)*(deriv_order+1);
if isnumeric(options.derivs) && options.derivs
    x=dde_select_tau(x,2,options.derivs,deriv_order);
end
nx=size(x,1);
nf=size(funcs.lhs_matrix(nx),1);
switch options.output
    case 'x'
        if nvec==0
            drhs=zeros(nf,nx,nvec);
            return
        end
        dx=reshape(full(sparse_blkdiag(repmat(eye(nx),1,1,ntaux))),nx,ntaux,nx*ntaux);
        dxvec=repmat(dx,1,1,nvec);
        dpvec=zeros(1,npar,nx*ntaux*nvec);
        xx=reshape(permute(x(:,:,:,ones(1,nx*ntaux)),[1,2,4,3]),nx,ntaux,nx*ntaux*nvec);
        pm=repmat(p,1,1,nx*ntaux);
        drhs=reshape(funcs.sys_dirderi{1}(xx,pm,dxvec,dpvec),[],nx,ntaux,nvec);
    case {'p','par','parameter'}
        nfp=length(options.free_par);
        if nfp==0
            drhs=zeros(nf,0,nvec);
            return
        end
        dpvec = dp_eye(npar,options.free_par,nvec);
        dxvec=zeros(nx,ntaux,nfp*nvec);
        xx=reshape(permute(x(:,:,:,ones(1,nfp)),[1,2,4,3]),nx,ntaux,nfp*nvec);
        pm=repmat(p,1,1,nfp);
        drhs=reshape(funcs.sys_dirderi{1}(xx,pm,dxvec,dpvec),[],nfp,nvec);
end        
end
function devp = dp_eye(npar,free_par, nvec)
nfp=length(free_par);
devp=zeros(npar,nfp);
devp(free_par,:)=eye(nfp);
devp=repmat(reshape(devp,1,npar,nfp),1,1,1,nvec);
devp=reshape(devp,1,npar,[]);%(1,npa,nfp*nv)
end