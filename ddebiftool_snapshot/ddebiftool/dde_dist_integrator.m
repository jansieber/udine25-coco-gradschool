function msh=dde_dist_integrator(cdeg,nint,type)
if nargin<3
    type='cheb';
end
switch type
    case {'cheb','finite','chebychev'}
        if length(nint)==1 && nint==round(nint)&&nint>0
            tcoarse=linspace(0,1,nint+1);
        else
            tcoarse=nint;
        end
        t_grid=dde_coll_meshfill(tcoarse,cdeg,'purpose','storage','grid','cheb','acc',true);
        is_cheb=true;
    case {'laguerre','infinite'}
        [t_grid,wt]=lagpts(cdeg,nint);
        t_grid=[0,t_grid(:).'];
        wt=[0,wt(:).'];
        is_cheb=false;
end
Jd=dde_coll_eva(zeros(0,length(t_grid)),t_grid,t_grid(2:end),cdeg,'diff',1,'output','matrix','submesh_limit',1);
Jd(end+1,1)=1;
Jd=Jd([end,1:end-1],:);
Ji=sparse(inv(full(Jd)));
if is_cheb
    wt=[0,full(Ji(end,2:end))];
end
ntau=length(t_grid);
tau=1:ntau-1;
msh=struct('Jibase',Ji,'t',t_grid,'wt',wt,'nint',nint,'degree',cdeg,...
    'laguerre',~is_cheb,'laguerre_alpha',nint,...
    'type',type,'ntau',ntau,'tau',tau);
msh.wJ={};
end