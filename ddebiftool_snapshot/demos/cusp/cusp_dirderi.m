function y = cusp_dirderi(order,xx,par,dx,dpar)
%% first and 2nd order directional derivatives in xx and par for cusp demo (see cusp_demo.m for equations)
%
% $Id$
%
%%
u=xx(1,2,:);
alpha=1/(1+exp(-4*u))-1/2;
dalpha=(4*exp(-4*u))./(exp(-4*u) + 1).^2;
ddalpha=(32*exp(-8*u))./(exp(-4*u) + 1).^3 -...
    (16*exp(-4*u))./(exp(-4*u) + 1).^2;
switch order
    case 1
        y=[-dx(1,1,:)+par(1)*dalpha.*dx(1,2,:)+dpar(1)*alpha-...
            par(2)*dx(2,2,:)-dpar(2)*xx(2,2,:)+dpar(4);...
            -dx(2,1,:)+par(3)*dalpha.*dx(1,2,:)+dpar(3)*alpha+dpar(5)];
    case 2
        y=[par(1)*ddalpha.*dx(1,2,:).^2+2*dpar(1)*dalpha.*dx(1,2,:)-...
            2*dpar(2)*dx(2,2,:);...
            par(3)*ddalpha.*dx(1,2,:).^2+2*dpar(3)*dalpha.*dx(1,2,:)];
end
end