function sys_tau_out=dde_tauin_from_tauseq(sys_tau_in,doreverse)
if nargin<2
    sys_tau_out=cellfun(@(c,i){i*ones(1,length(c))},...
        sys_tau_in,num2cell(1:length(sys_tau_in)));
    if ~iscell(sys_tau_out)
        sys_tau_out={sys_tau_out};
    end
    sys_tau_out=cat(2,sys_tau_out{:});
elseif ismember(doreverse,{'doreverse','reverse'})
    ts=find(diff([0,sys_tau_in]));
    te=find(diff([sys_tau_in,sys_tau_in(end)+1]));
    sys_tau_out=arrayfun(@(s,e){s:e},ts,te);
end
end
