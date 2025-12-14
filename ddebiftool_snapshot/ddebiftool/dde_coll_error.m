function [out1,out2,out3]=dde_coll_error(coll,varargin)
%% interpolation error estimate for piecewise polynomials
% based on difference between polynomial interpolations across subinterval
% boundaries (so, for single polynomials this is always 0)
default={'output','solution'};
options=dde_set_options(default,varargin,'pass_on');
err=coll;
err.profile=NaN(size(err.profile));
if size(coll.mesh,2)<=coll.degree+1
    [out1,out2,out3]=outputs(options.output,err);
    return
end
x=coll.profile;
t=coll.mesh;
d=coll.degree;
t_int=cat(2,1.5*t(1)-0.5*t(2),0.5*(t(1:end-1)+t(2:end)),1.5*t(end)-0.5*t(end-1));
x_int=dde_coll_eva(x,t,t_int,d,'assert_boundaries',false);
err.profile(:,1:d)=dde_coll_eva(x_int(:,1:d+2),t_int(1:d+2),t(1:d),d+1);
err.profile(:,end-d+1:end)=dde_coll_eva(x_int(:,end-d-1:end),t_int(end-d-1:end),t(end-d+1:end),...
    d+1);
err.profile(:,d+1:end-d)=dde_coll_eva(x_int(:,1:end-1),t_int(1:end-1),t(d+1:end-d),d);
err.profile=err.profile-x;
[out1,out2,out3]=outputs(options.output,err);
end
function [out1,out2,out3]=outputs(output,err)
d=err.degree;
err_comp=max(abs(err.profile),[],1);
err_max=max(err_comp);
err_comp=cat(1,reshape(err_comp(1:end-1),d,[]),err_comp(d+1:d:end));
err_comp=max(err_comp,[],1);
meas=cumsum([0,err_comp.^(1/(d+1))]);
switch output
    case {'max','maximum','norm'}
        [out1,out2,out3]=deal(err_max,err,meas);
    case {'measure'}
        [out1,out2,out3]=deal(meas,err,err_max);
    otherwise
        [out1,out2,out3]=deal(err,err_max,meas);
end        
end