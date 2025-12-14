function pout=dde_point_embed(points,template_inp,free_par_ind,xind)
template=feval(['dde_',template_inp(1).kind,'_create'],template_inp(1));
template=dde_trim_point(template,template_inp(1));
pout=repmat(template,size(points));
if isempty(points)
    return
end
if nargin<3
    free_par_ind=1:size(template.parameter,2);
end
[fields,xfields]=feval(['dde_',template.kind,'_vars']);
if nargin<4
    xind=1:size(template.(xfields{1}),1);
end
npnames=setdiff(fieldnames(fields),'parameter');
nxnames=setdiff(npnames,xfields);
for i=1:npts
    pout(i).parameter(free_par_ind)=points(i).parameter(1:length(free_par_ind));
    for k=1:length(nxnames)
        pout(i).(nxnames{k})=points(i).(nxnames{k});
    end
    for k=1:length(xfields)
        template_xk=template.(xfields{k});
        template_xk(xind,:)=points(i).(xfields{k});
        pout(i).(xfields{k})=template_xk;
    end
end
end
