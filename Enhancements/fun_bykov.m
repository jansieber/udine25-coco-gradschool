function y=fun_bykov(ord,iu,ip,u,p,du,dp)
[    Q1,        Q2,        Q3,        Q4,        Q5,        Q6,        K]=deal(...
p(ip.Q1,:),p(ip.Q2,:),p(ip.Q3,:),p(ip.Q4,:),p(ip.Q5,:),p(ip.Q6,:),p(ip.K,:));
[      x,        y,        s]=deal(...
  u(iu.x,:),u(iu.y,:),u(iu.s,:));
z  = 1-x-y-s;
if ord==0
    xp = 2*Q1.*z.^2-2*Q5.*x.^2-Q3.*x.*y;
    yp = Q2.*z-Q6.*y-Q3.*x.*y;
    sp = Q4.*z-K.*Q4.*s;
    y=[xp;yp;sp];
    return
end
[    dQ1,        dQ2,        dQ3,        dQ4,        dQ5,        dQ6,        dK]=deal(...
dp(ip.Q1,:),dp(ip.Q2,:),dp(ip.Q3,:),dp(ip.Q4,:),dp(ip.Q5,:),dp(ip.Q6,:),dp(ip.K,:));
[      dx,        dy,        ds]=deal(...
  du(iu.x,:),du(iu.y,:),du(iu.s,:));
dz  = -dx-dy-ds;
dxp=2*dQ1.*z.^2+4*Q1.*z.*dz-2*dQ5.*x.^2-4*Q5.*x.*dx...
    -dQ3.*x.*y-Q3.*dx.*y-Q3.*x.*dy;
dyp=dQ2.*z+Q2.*dz-dQ6.*y-Q6.*dy-dQ3.*x.*y-Q3.*dx.*y-Q3.*x.*dy;
dsp=dQ4.*z+Q4.*dz-dK.*Q4.*s-K.*dQ4.*s-K.*Q4.*ds;
y=[dxp;dyp;dsp];
end
