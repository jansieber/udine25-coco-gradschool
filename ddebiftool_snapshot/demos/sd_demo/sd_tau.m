%% User-provided state-dependent delays
% Implementation for tutorial sd_demo:
%
%   function dtau=sd_tau(k,xx,par)
%
% returns delay $\tau_k$, depending on |xx(:,1:k)|
% ($[x(t),\ldots,x(t-\tau_{k-1})]$) and |par|.
%
% As tau2 does not depend on x(t-tau1) tau1, tau2 can be returned at the
% same time, controlled by sys_tau_seq. Similarly, tau3,...tau6 can all be
% returned together as they do not dependon x(t-tau3),...
%%
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
% if k==1
%     tau=p(10)+pad;
% elseif k==2
%     tau=p(11)+pad;
% elseif k==3
%     tau=2+p(5)*p(10)*x(2,1,:).*x(2,2,:);
% elseif k==4
%     tau=1-1./(1+x(2,3,:).*x(1,1,:));
% elseif k==5
%     tau=x(4,1,:);
% elseif k==6
%     tau=x(5,1,:);
% end
% end
% 
