function [in,notin]=dde_options_filter(prefix,args,symbol)
if nargin<3
    symbol='.';
end
args=args(:)';
names=args(1:2:end-1);
values=args(2:2:end);
nf={'UniformOutput',false};
if isempty(prefix)
    symcount=strfind(anames,symbol);
    hassym=~cellfun(@isempty,symcount);
    innames=cellfun(@(s,i)s(i{1}:end),names(hassym),symcount(hassym),nf{:});
else
    ps=[prefix,symbol];
    lenps=length(ps);
    hassym=strncmp(names,ps,lenps);
    innames=cellfun(@(s)s(lenps+1:end),names(hassym),nf{:});
end
invals=values(hassym);
in=reshape(cat(1,innames,invals),1,[]);
notin=reshape(cat(1,names(~hassym),values(~hassym)),1,[]);
end
