function [prob,modids] = coll_construct_modes(prob, collfid,varargin)
default={'orders',1,'mid','mode','norm',true,'fid','modes'};
options=sco_set_options(default,varargin,'pass_on');
efid = coco_get_id(collfid, options.fid);
[colldata,uidx]=coco_get_func_data(prob,collfid,'data','uidx');
seg  = colldata.coll_seg;
[data,modids]=init_data(collfid,seg,options.orders,options.mid,options.norm);
prob = coco_add_func(prob, efid, @coll_modes, data, ...
  'inactive', modids, 'uidx', uidx(seg.maps.xbp_idx), ...
  'remesh', @coll_modes_remesh,'F+DF');
end
function [data,modids]=init_data(collfid,seg,orders,basename,incl_norm)
xbp_shp=seg.maps.xbp_shp;
dim=xbp_shp(1);
no=length(orders);
[nmo,nmdim]=meshgrid(orders,1:dim);
names=arrayfun(@(o,n)sprintf('%s_%d_%d',basename,o,n),nmo,nmdim,...
    'uniformoutput',false);
if incl_norm
    names=[{[basename,'_sum']},names(:)'];
end
modids=coco_get_id(collfid,names);
factors=sqrt(2)*(orders~=0)+1*(orders==0);
factors=reshape(repmat(factors(:)',2,1),[],1);
facmat=diag(factors);
data=struct('fid',collfid,'orders',orders,'facmat',facmat,'names',{names},...
    'incl_norm',incl_norm);
end
function [data,y,J] = coll_modes(prob, data, u)
colldata=coco_get_func_data(prob,data.fid,'data');
seg = colldata.pr.coll_seg;
maps = seg.maps;
msh  = seg.mesh;
xbp_shp=seg.maps.xbp_shp;
dim=xbp_shp(1);
x     = u;
no=length(data.orders);
frmat=cat(3,cos(2*pi*msh.tbp*data.orders),sin(2*pi*msh.tbp*data.orders));
frmat=reshape(permute(frmat,[1,3,2]),[],2*no);
mat=kron(data.facmat*frmat.',eye(dim));
mat=reshape(mat,dim,2,no,[]);
mat=reshape(permute(mat,[2,1,3,4]),dim*2*no,[]);
intJ=mat*maps.W'*msh.wts2*msh.kas2*maps.W*(0.5)/maps.NTST;
smat=kron(eye(no*dim),[1,1]);
Jx=intJ*x;
yf=smat*(Jx).^2;
Jf=smat*2*diag(Jx)*intJ;
ya=NaN(0,1);
Ja=NaN(0,size(Jf,2));
if data.incl_norm
    ya=sqrt(sum(yf));
    Ja=sum(Jf,1)/(2*ya);
end 
y=[ya;yf];
J=[Ja;Jf];
end
function [prob, stat, xtr] = coll_modes_remesh(prob, data, chart, ub, Vb) %#ok<INUSD>
colldata=coco_get_func_data(prob,data.fid,'data');
seg  = colldata.coll_seg;
maps = seg.maps;
xtr  = [];
uidx = coco_get_func_data(prob, seg.fid, 'uidx');
prob = coco_change_func(prob, data, 'uidx', uidx(maps.xbp_idx));
stat = 'success';
end
