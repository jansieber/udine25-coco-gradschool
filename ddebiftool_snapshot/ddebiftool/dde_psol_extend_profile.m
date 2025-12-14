function [tmsh_ext,prof_ext]=dde_psol_extend_profile(pt,t,varargin)
%% extend periodic solution profile and mesh beyond the interval tmesh
% for stability computations
% extension is by multiples of base collocation intervals
% extension of the profile (optional) permits difference between
% profile(:,1) and profile(:,end). This difference will be added to each
% extension.
default={{'c','collocation_parameters'},[]};
options=dde_set_options(default,varargin,'pass_on');
[tmsh,deg]=deal(pt.mesh,pt.degree);
bdmsh=tmsh([1,end]);       % mesh interval boundaries
lenmsh=tmsh(end)-tmsh(1);  % mesh interval length
mrect=cat(1,reshape(tmsh(1:end-1),deg,[]),reshape(tmsh(deg+1:deg:end),1,[]));
tmin=min([t(:);bdmsh(1)]);
tmax=max([t(:);bdmsh(2)]);
%% hack to provide history for smoothness condition if force_smooth
% is enabled but delays are all 0. This should add only zero Floquet
% multipliers.
if tmin==bdmsh(1) && strcmp(options.c,'force_smooth')
    tmin=tmin+mrect(1,end)-mrect(end,end);
end
npastmeshes=ceil((bdmsh(1)-tmin)/lenmsh);
nfutmeshes=ceil((tmax-bdmsh(2))/lenmsh);
nmeshes=npastmeshes+1+nfutmeshes;
mrect_rep=reshape(mrect(:,:,ones(nmeshes,1)),deg+1,[]);
mrect_add=reshape(lenmsh*(-npastmeshes:nfutmeshes),1,1,nmeshes);
mrect_add=reshape(mrect_add(ones(deg+1,1),ones(size(mrect,2),1),:),deg+1,[]);
mrect_ext=mrect_rep+mrect_add;
m_ext=cat(2,reshape(mrect_ext(1:end-1,:),1,[]),mrect_ext(end,end));
itmin=find(mrect_ext(1,:)<=tmin,1,'last');
itmax=find(mrect_ext(end,:)>=tmax,1,'first');
mrect_cut=mrect_ext(:,itmin:itmax);
tmsh_ext=[reshape(mrect_cut(1:end-1,:),1,[]),mrect_cut(end,end)];
if nargout==1
    return
end
itmin_ext=find(m_ext==tmsh_ext(1));
itmax_ext=find(m_ext==tmsh_ext(end));
profile=pt.profile;
n=size(profile,1);
profile_diff=profile(:,end)-profile(:,1);
profile_rep=profile(:,1:end-1,ones(nmeshes,1));
profile_add=reshape(profile_diff*(-npastmeshes:nfutmeshes),n,1,nmeshes);
profile_add=profile_add(:,ones(1,size(profile_rep,2)),:);
profile_all=profile_rep+profile_add;
profile_final=profile_all(:,1,end)+profile_diff;
profile_all=[reshape(profile_all,n,[]),profile_final];
prof_ext=profile_all(:,itmin_ext:itmax_ext);
end
