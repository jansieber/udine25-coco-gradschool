function [points,ind]=dde_point_from_x(x,template_inp,free_par_ind,xind)
%% insert variable values into point structure
%
% $Id: dde_point_from_x.m 308 2018-10-28 15:08:12Z jansieber $
%%
template=feval(['dde_',template_inp(1).kind,'_create'],template_inp(1));
template=dde_trim_point(template,template_inp(1));
if nargin<3
    free_par_ind=1:size(template.parameter,2);
end
[ind,len,xfields]=dde_ind_from_point(template,free_par_ind); %#ok<ASGLU>
if nargin<4
    xind=1:size(template.(xfields{1}),1);
end
sx=size(x);
ptshape=[sx(2:end),1];
points=repmat(template,ptshape);
if isempty(x)
    return
end
npnames=setdiff(fieldnames(ind),'parameter');
nxnames=setdiff(npnames,xfields);
npts=prod(ptshape);
x=reshape(x,size(x,1),npts);
for i=1:npts
    points(i).parameter(free_par_ind)=x(ind.parameter,i);
    for k=1:length(nxnames)
        fd=ind.(nxnames{k});
        if ~isstruct(fd)
            points(i).(nxnames{k})=reshape(x(fd,i),size(template.(nxnames{k})));
        else
            y=x(reshape(fd.re,[],1),i)+1i*x(reshape(fd.im,[],1),i);
            points(i).(nxnames{k})=reshape(y,size(template.(nxnames{k})));
        end        
    end
    for k=1:length(xfields)
        fd=ind.(xfields{k});
        pfield=points(i).(xfields{k});
        template_xk=template.(xfields{k});
        if ~isstruct(fd)
            pfield(xind,:)=reshape(x(fd,i),size(template_xk(xind,:)));
        else
            y=x(reshape(fd.re,[],1),i)+1i*x(reshape(fd.im,[],1),i);
            pfield(xind,:)=reshape(y,size(template_xk(xind,:)));
        end
        points(i).(xfields{k})=pfield;
    end
end
end
