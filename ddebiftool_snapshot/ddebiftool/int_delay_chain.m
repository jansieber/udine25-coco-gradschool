function [res,yarr,ind]=int_delay_chain(ind,order,xx,p,bd,x0,t,val,dx,dp,dbd,dx0,dt,dval)
%% approximate algebraic equation for distributed delay
% composed of sequence nested integrals
%
% y_k(t-r)=x0(t)+int_0^r g_k(s,x_k(t-s),p) dr for r=0..tau
%
% (x0(t) can be dynamic variable or parameter) discretized
%
% y_k=intwJ*[x0;tau*g(t(2:end)*tau,xx_k(:,t-t(2:end)*tau),p)]
%
% where wt(i) and t(i) are predefined integration weights and nodes on the
% interval [0,1], if the interval is finite, or on [0,infty] if the
% interval is infinite (Laguerre or generalised Laguerre approximation)
%
% xx_k(j,:) can be defined as y_i(l,:) for some i<k and l<=size(y_i,1)
% through the index list ind.xinputs
igx=ind.args.igx.local.x{1};
yarr=int_delay_chain_link(ind.histories,ind.args.igp,ind.msh,order,...
    xx(igx,:,:),p,bd,x0, dx(igx,:,:),dp,dbd,dx0);
res=int_delay_chain_sum(ind.sum,ind.msh,order,yarr,bd,val,t,dbd,dval,dt);
end
%%
function yarr=int_delay_chain_link(hist,igp,msh,order,varargin)
nbaseargs=4;
args=cell(nbaseargs,length(varargin)/nbaseargs);
args(1:length(varargin))=varargin;
[xx,p,bd,x0,gargc]=deal(cell(1,order+1));
for i=1:min(2,order+1)
    [xx{i},p{i},bd{i},x0{i}]=deal(args{:,i});
end
msh.t=msh.t;
msh.int=hist.int;
msh.id=hist.id;
ordlist=0:order;
[nw,nvec,nord]=deal(length(msh.t),size(bd{1},2),length(ordlist));
yarr=NaN(hist.lastxarg,nw,nvec,nord);
yarr(hist.xvals{1}{hist.id},:,:,:)=0;
for i=1:min(2,order+1)
    yarr(hist.xvals{1}{hist.id},:,:,i)=...
        reshape(xx{i},length(hist.xvals{1}{hist.id}),nw,nvec);
end
for k=1:length(hist.fcn)
    for i=1:order+1
        gargc{i}=int_delay_chain_ind2args(hist.xinputs,hist.xinit,...
            igp,i-1,k,  yarr(:,:,:,i),p{i},bd{i},x0{i});
    end
    msloc=setfield(msh,'wJ',msh.wJ{k}); %#ok<SFLD>
    gargs=[gargc{:}];
    out=int_delay_chain_element(hist.fcn{k},msloc,ordlist,gargs{:});
    int=hist.xvals{k+1}{hist.int};
    yarr(int,:,:,:)=cat(4,out{hist.int,:});
    id=hist.xvals{k+1}{hist.id};
    yarr(id,:,:,:)=cat(4,out{hist.id,:});
end
end
%%
function args=int_delay_chain_ind2args(ixc,ixinitc,igp,order,ifun,xvals,p,bd,x0)
%% construct arguments for delay integrands from ind structure
ixvals=ixc{ifun};
ixinit=ixinitc{ifun};
args{1}=xvals(ixvals,:,:);
if order<2
    args(2:4)={p(1,igp.local.parameter{ifun},:),bd,x0(ixinit,:)};
end
end
