function [rts,rtval,rtprime]=dde_coll_roots(pt,c,varargin)
%% find roots of c*x(t)-res or c*x'(t)-res on interval [0,1]
% (may fail for phase rotations, as it assumes periodicity)
default={'diff',0,'safety',1e-3,'res',0,'relative',true,'vectorized',false};
options=dde_set_options(default,varargin,'pass_on');
%% extend coll periodically by one subinterval (accommodating rotations)
pt_orig=pt;
pt.mesh=[pt.mesh,pt.mesh(end)+pt.mesh(2:pt.degree+1)];
ondeg=ones(1,pt.degree);
pt.profile=[pt.profile,pt.profile(:,end)*ondeg-pt.profile(:,1)*ondeg+pt.profile(:,2:pt.degree+1)];
%% find subintervals and values of c*x-res or c(x,p)
ic=1:pt.degree:length(pt.mesh);
nint=length(ic)-1;
if isnumeric(c)
    pt.profile=c*pt.profile-options.res;
elseif isa(c,'function_handle') && options.vectorized
    pt.profile=c(pt.profile,pt.parameter)-options.res;
elseif isa(c,'function_handle') && ~options.vectorized
    for i=size(pt.profile,2):-1:1
        profile(i)=c(pt.profile(:,i),pt.parameter)-options.res;
    end
    pt.profile=profile;
end
if options.diff>0
    d_opts={'diff',options.diff,'kron',true};
    cxvals0=dde_coll_eva(pt,pt.mesh(1:end-1),d_opts{:})-options.res;
    cxvals1=dde_coll_eva(pt,pt.mesh(ic(2:end)),d_opts{:},'submesh_limit',1)-options.res;
else
    cxvals_all=pt.profile;
    cxvals1=cxvals_all(pt.degree+1:pt.degree:end);
    cxvals0=cxvals_all(1:end-1);
end
%% each column is a subintervals and the corresponding time points
tvals=cat(1,reshape(pt.mesh(1:end-1),pt.degree,nint),pt.mesh(ic(2:end)));
cxvals=cat(1,reshape(cxvals0,pt.degree,nint),cxvals1);
%% find which intervals have sign changes
maxcxvals=max(cxvals,[],1);
mincxvals=min(cxvals,[],1);
icross=sign(maxcxvals.*mincxvals)<0;
%% for (potentially discontinous) derivative find sign changes between intervals
ijump=[sign(cxvals(end,1:end-1).*cxvals(1,2:end))<0,false];
%% check also for intervals with small values
if options.relative
    margin=max(abs(cxvals(:)))*options.safety;
else
    margin=options.safety;
end
ismall=any(abs(cxvals)<margin,1);
icheck=find(icross|ismall);
rtc=cell(1,length(icheck)+1);
rtc{end}=tvals(end,ijump);
for i=1:length(icheck)
    rtc{i}=dde_coll_local_roots(tvals(:,icheck(i))',cxvals(:,icheck(i)).',2);
    rtc{i}=rtc{i}(:).';
end
rts=[rtc{:}];
rts=rts(rts<=pt_orig.mesh(end));
rtval=dde_coll_eva(pt,rts);
rtprime=dde_coll_eva(pt,rts,'diff',1);
end
function r=dde_coll_local_roots(tvals,xvals,restrict)
%% Find roots of polynomial stored at tvals with values xvals 
% uses chebyshev formula for polynomial root finding based on eig third
% argument has values 0,1,2 (default 2) restrict>0 restricts to real roots,
% restrict>1 restricts to roots between tvals(1) and tvals(end) tvals is
% assumed to be ordered (increasing or decreasing) and
% tvals(end)=\=tvals(1).
t=2*(tvals-tvals(1))/(tvals(end)-tvals(1))-1;
d=length(t)-1;
t=t(:)';
Aint=[ones(1,d+1);t;NaN(d-1,d+1)];
for i=3:d+1
    Aint(i,:)=2*t.*Aint(i-1,:)-Aint(i-2,:);
end
ccheb=xvals/Aint;
if all(abs(ccheb)==0)
    iclead=2;
else
    iclead=find(abs(ccheb)~=0,1,'last');
end
c=ccheb(:,1:iclead-1)./repmat(ccheb(iclead),1,iclead-1);
d=iclead-1;
cm0=diag([ones(d-2,1);2],-1)+diag(ones(d-1,1),1);
nx=size(xvals,1);
r=cell(nx,1);
if nargin<3
    restrict=2;
end
for i=1:nx
    cm=cm0;
    cm(1,:)=cm(1,:)-c(end:-1:1);
    r{i}=eig(cm/2);
    if restrict>0
        r{i}=r{i}(imag(r{i})==0);
    end
    if restrict>1
        r{i}=r{i}(r{i}<=1&r{i}>=-1);
    end
    r{i}=(r{i}+1)*0.5*(tvals(end)-tvals(1))+tvals(1);
end
if nx==1
    r=r{1};
end
end