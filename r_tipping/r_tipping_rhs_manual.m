function f=r_tipping_rhs_manual(order,ip,y,pa,dy,dpa)
%% r.h.s. for R-tipping problem
% assign names to parameters and variables
[      r,        cm,         beta,         omega,         cp,         a]=deal(...
 pa(ip.r,:),pa(ip.cm,:),pa(ip.beta,:),pa(ip.omega,:),pa(ip.cp,:),pa(ip.a,:)); % parameter names
[x,p,q]=deal(y(1,:),y(2,:),y(3,:));
%% r.h.s
rho=p.^2+q.^2;
srho=sqrt(rho);
L=cm+rho.*(cp-cm)+beta.*q.*srho;
nlin=r-r.*rho;
if order==0
    f=[(x+L).^2-a;   p.*nlin+omega.*q;  -omega.*p+q.*nlin];
    return
end
%% assign names to deviations
[     dr,        dcm,         dbeta,         domega,         dcp,         da]=deal(...
dpa(ip.r,:),dpa(ip.cm,:),dpa(ip.beta,:),dpa(ip.omega,:),dpa(ip.cp,:),dpa(ip.a,:)); % parameter names
[dx,dp,dq]=deal(dy(1,:),dy(2,:),dy(3,:));
%% 1st order directional derivative
drho=2*(dp.*p+dq.*q);
dsrho=drho./srho/2;
dL=dcm+rho.*(dcp-dcm)+drho.*(cp-cm)+dbeta.*q.*srho+beta.*dq.*srho+beta.*q.*dsrho;
dnlin=dr-dr.*rho-r.*drho;
if order==1
    f=[2*(x+L).*(dx+dL)-da; dp.*nlin+p.*dnlin+domega.*q+omega.*dq;...
        -domega.*p-omega.*dp+dq.*nlin+q.*dnlin];
    return
end
%% 2nd order directional derivative
d2rho=2*(dp.^2+dq.^2);
d2srho=(dp.*q-dq.*p).^2./srho.^3;
d2L=d2rho.*(cp-cm)+2*drho.*(dcp-dcm)+2*dbeta.*dq.*srho+2*beta.*dq.*dsrho+2*dbeta.*q.*dsrho+...
    beta.*q.*d2srho;
d2nlin=-2*dr.*drho-r.*d2rho;
f=[2*(dx+dL).^2+2*(x+L).*d2L; 2*dp.*dnlin+p.*d2nlin+2*domega.*dq;...
    -2*domega.*dp+2*dq.*dnlin+q.*d2nlin];
end
