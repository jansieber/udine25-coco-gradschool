function pout=dde_point_restrict(points,free_par_ind,xind)
pout=points;
if isempty(points)
    return
end
if nargin<3
    free_par_ind=1:size(points(1).parameter,2);
end
[fields,xfields]=feval(['dde_',points(1).kind,'_vars']); %#ok<ASGLU>
if nargin<4
    xind=1:size(points(1).(xfields{1}),1);
end
for i=1:npts
    pout(i).parameter=points(i).parameter(free_par_ind);
    for k=1:length(xfields)
        template_xk=points(i).(xfields{k});
        template_xk=template_xk(xind,:);
        pout(i).(xfields{k})=template_xk;
    end
end
end
