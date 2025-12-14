function args=int_delay_chain_ind2args(ixinputs,ipinputs,order,ifun,xargs,p,bd)
%% construct arguments for delay integrands from ind structure
ix=ixinputs{ifun};
args{1}=xargs(ix,:,:);
if order<2
    args{2:3}={p(1,ipinputs.map.parameter{ifun},:),bd};
end
end