function yarr=int_delay_chain_link(ind,order,varargin)
nbaseargs=3;
args=cell(nbaseargs,order+1);
args(1:length(varargin))=varargin;
[xx,p,bd,gargc]=deal(cell(1,ord+1));
for i=1:ord+1
    [xx{i},p{i},bd{i}]=deal(args{:,i});
end
msh.t=ind.msh.t;
msh.intid=ind.intid;
ordlist=0:order;
[nw,nvec,nord]=deal(length(msh.t),size(bd{1},2),length(ordlist));
yarr=NaN(ind.lastxarg,nw,nvec,nord);
for k=1:length(ind.fcn)
    for i=1:order+1
        gargc{i}=int_delay_chain_ind2args(i-1,k,ind,yarr(:,:,:,i),p{i},bd{i});
    end
    msh.wJ=ind.msh.wJ{k};
    gargs=[gargc{:}];
    out=int_delay_chain_element(ind.fcn{k},msh,ordlist,gargs{:});
    int=ind.xargs{k+1}{ind.intid.int};
    yarr(int,:,:,:)=cat(4,out{:,ind.intid.int});
    id=ind.xargs{k+1}{ind.intid.id};
    yarr(id,:,:,:)=cat(4,out{:,ind.intid.id});
end
end
%%
function args=int_delay_chain_ind2args(order,ifun,ind,xargs,p,bd)
%% construct arguments for delay integrands from ind structure
ix=ind.xinputs{ifun};
args{1}=xargs(ix,:,:);
if order<2
    args{2:3}={p(1,ind.pinputs{ifun},:),bd};
end
end
