function [points,farey]=dde_locate_resonances(trfuncs,trbranch,n)
%% determine all resonanes up to order n along branch of torus bifurcations
%% determine pseudo-arclength along branch
np=length(trbranch.point);
arclength=zeros(1,np);
for i=2:np
    pdiff=p_axpy(-1,trbranch.point(i-1),trbranch.point(i));
    arclength(i)=p_norm(pdiff);
end
arclength=cumsum(arclength);
%% arclength is used as base parameter to interpolate points on branch with splines
rotation=trfuncs.get_comp(trbranch.point,'omega')/2;
maxrot=ceil(max(rotation));
minrot=floor(min(rotation));
basefarey=dde_farey(n);
allfarey=basefarey;
for i=2:maxrot
    allfarey=[allfarey,[basefarey(1,2:end-1)+(i-1)*basefarey(2,2:end-1);...
        basefarey(2,2:end-1)]];
end
for i=-1:-1:minrot
    allfarey=[[-basefarey(1,end:-1:2)+(i+1)*basefarey(2,end:-1:2);...
        basefarey(2,end:-1:2)],allfarey];
end
fareyfrac=allfarey(1,:)./allfarey(2,:);
% omega is between 0 and 1, 1 corresponds to 1:2 resonance
% for i=1:length(rotation)
%     if rotation(i)<0 ||rotation(i)>1
%         extrarot=floor(rotation(i));
%         ev=trfuncs.get_comp(trbranch.point(i),'eigenvector');
%         dim=size(ev.profile,1)/2;
%         evz=ev.profile(1:dim,:)+1i*ev.profile(1:dim,:);
%         fac=exp(1i*2*pi*extrarot*ev.mesh);
%         evz=evz.*repmat(fac,[dim,1]);
%         trbranch.point(i).profile(dim+(1:2*dim),:)=[real(evz);imag(evz)];
%         rotation(i)=mod(rotation(i),1);
%         trbranch.point(i).parameter(end-1)=rotation(i)*2;
%     end
% end
fareyselect=fareyfrac<=max(rotation)&fareyfrac>=min(rotation);
allfarey=allfarey(:,fareyselect);
fareyfrac=fareyfrac(fareyselect);
%% interpolate to find rotation for fareyfrac on branch
resonances=[];
resonance_values=[];
farey=[];
for i=2:np
    select=find(sign((fareyfrac-rotation(i-1)).*(fareyfrac-rotation(i)))<=0);
    fsel=fareyfrac(select);
    if isempty(fsel)
        continue
    end
    if rotation(i)==rotation(i-1)
        resonances=[resonances,arclength(i)]; %#ok<AGROW>
        resonance_values=[resonance_values,fsel]; %#ok<AGROW>
    else
        slope=(arclength(i)-arclength(i-1))/(rotation(i)-rotation(i-1));
        arc=arclength(i-1)+(fsel-rotation(i-1))*slope;
        if length(fsel)>1 && arc(end)<arc(1)
            fsel=fsel(end:-1:1);
            select=select(end:-1:1);
            arc=arc(end:-1:1);
        end
        resonances=[resonances,arc]; %#ok<AGROW>
        resonance_values=[resonance_values,fsel]; %#ok<AGROW>
        farey=[farey,allfarey(:,select)]; %#ok<AGROW>
    end
end
nres=length(resonances);
%% return interpolated points
points=repmat(trbranch.point(1),[1,nres]);
free_par_ind=trbranch.parameter.free;
[dim,nt]=size(trbranch.point(1).profile);
extract=@(f)cell2mat(arrayfun(f,trbranch.point,'uniformoutput',false));
pars=extract(@(x)x.parameter(free_par_ind)')';
periods=[trbranch.point.period];
meshes=extract(@(x)x.mesh')';
profiles=reshape(extract(@(x)x.profile),[dim,nt,np]);
res_pars=interp1(arclength,pars,resonances,'spline');
res_periods=interp1(arclength,periods,resonances,'spline');
res_meshes=interp1(arclength,meshes,resonances,'spline');
res_profiles=interp1(arclength,permute(profiles,[3,1,2]),resonances,'spline');
for i=1:nres
    points(i).parameter(free_par_ind)=res_pars(i,:);
    points(i).parameter(free_par_ind(end-1))=resonance_values(i)*2;
    points(i).period=res_periods(i);
    points(i).mesh=res_meshes(i,:);
    points(i).profile=reshape(res_profiles(i,:,:),dim,nt);
end
end

