%% User-provided state-dependent delays for tutorial sd_demo:
function tau=sd_tau(k,x,p)
if all(k<3)
    tau=cat(1,p(1,10,:),p(1,11,:));
    tau=tau(k,:);
else
    tau=cat(1,...
        2+p(1,5,:).*p(1,10,:).*x(2,1,:).*x(2,2,:),...
        1-1./(1+x(2,3,:).*x(1,1,:)),...
        x(4,1,:),...
        x(5,1,:));
    tau=tau(k-2,:);
end
end
