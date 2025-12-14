function y=int_delay_chain_sum(gsum,msh,order,yarr,bd,val,t,dbd,dval,dt)
yarr=yarr(gsum.input,:,:,:);
[ny,nw,nvec,nord]=size(yarr);
valc={val,dval,zeros(size(val,1),nvec)};
y=-valc{order+1};
%% integral from 0 to bd if requested (t=NaN)
nt=size(t,1);
tyrg=1:ny*nt;
nyt=length(tyrg);
yt=zeros(nyt,nvec);
yval=permute(reshape(repmat(yarr,[nt,1,1,1]),[ny,nt,nw,nvec,nord]),...
    [3,1,2,4,5]);
t_is_whole_interval=isnan(t(:,1));
iwhole=find(t_is_whole_interval,1,'first');
if any(t_is_whole_interval)
    whrg=sub2ind([ny,nt],1:ny,ones(1,ny)*iwhole);
    tyrg(whrg)=[];
    yval=reshape(yval(:,1,:,:,:),nw,ny*nvec*(order+1));
    yval=reshape(msh.wt(:).'*yval,ny,nvec,order+1);
    fac=gsum.coeffs_whole.count{order+1};
    cf= gsum.coeffs_whole.coeffs{order+1};
    bc=cat(3,bd,dbd);
    for k=1:length(fac)
        yt(whrg,:)=yt(whrg,:)+...
            (yval(:,:,cf(k,1)).*bc(:,:,cf(k,2)))*fac(k);
    end
end
if ~isempty(tyrg)
    %% integral from 0 to t*bd for variable relative boundaries
    nt=nt-length(iwhole);
    yval=reshape(yval,nw,nyt*nvec,order+1);
    tc={speye(nt*nvec),diag(sparse(dt(:))),diag(sparse(dt(:).^2))};
    Jw=cell(1,order+1);
    for k=1:order+1
        Jw{k}=dde_coll_eva(zeros(0,nw),msh.t,t(:).',msh.degree,...
            'output','matrix','diff',k-1,'assert_boundaries',false);
        Jw{k}=reshape(repmat(tc{k}*Jw{k},1,ny).',nw,nyt*nvec);
    end
    mJy=@(Ja,ya)reshape(sum(Ja.*ya,1),nyt,nvec);
    fac=gsum.coeffs_t.count{order+1};
    cf= gsum.coeffs_t.coeffs{order+1};
    for k=1:length(fac)
        sel=cf(k,:);
        yt(tyrg,:)=yt(tyrg,:)+mJy(Jw{sel(1)},yval(:,:,sel(2)))*fac(k);
    end
end
yt=reshape(permute(reshape(yt,ny,nt,nvec),[2,1,3]),nt*ny,nvec);
y=y+gsum.M*yt;
end
