function out=dde_coll_bd_diff(pt,varargin)
default={'boundary',false,'order',1,'output','profile',...
    {'collocation_parameters','c'},[],...
    'mesh',pt.mesh,'range',@(t)t>0&t<1,'force_smooth',false};
options=dde_set_options(default,varargin,'pass_on');
tcoarse=options.mesh(1:pt.degree:end);
t=zeros(2,0);
force_smooth=options.force_smooth||(ischar(options.collocation_parameters) &&...
    strcmp(options.collocation_parameters,'force_smooth'));
if force_smooth
    if options.boundary
        t=[tcoarse(end);tcoarse(1)];
    else
        t=[tcoarse(options.range(tcoarse));...
           tcoarse(options.range(tcoarse))];
    end
end
args={'diff',options.order,'kron',true,'output',options.output};
dt1=dde_coll_eva(pt,t(1,:),'submesh_limit',1,args{:});
dt0=dde_coll_eva(pt,t(2,:),'submesh_limit',0,args{:});
out=dt1-dt0;
end
