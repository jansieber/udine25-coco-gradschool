function [opts,pass_on]=loc_set_options(default,varargin,pass_on)
opts=struct(default{:});
default=reshape(default,2,[]);
args=reshape(varargin,2,[]);
dnames=default(1,:);
argnames=args(1,:);
[~,ism]=ismember(argnames,dnames);
assert(nargin>2 || ~all(ism),'loc_set_options: unknown option')
for i=1:length(ism)
    opts.(dnames{ism(i)})=args{2,i};
end
pass_on=reshape(args(:,~ism),1,[]);
end
