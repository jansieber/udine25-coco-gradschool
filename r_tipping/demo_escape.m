%% Demo showing phase dependence of escape for rate-induced tipping
clear
format compact
pnames={'r','cm','beta','omega','cp','a','phi'}; % parameter names
vnames={'x','p','q'};
[ip,npars]=structind_from_names(pnames);
[iv,dim]=structind_from_names(vnames);
[cm, beta, omega, cp,   a]=deal(...
  0, 0.5,   2*pi,  4, 1.75);
L=@(t,r,phi)Lfun(t,cm,beta,omega,cp,r,phi);
r_rg=2.1:0.05:2.8;
phi_rg=linspace(0,2*pi,10);
phi_rg=phi_rg(1:end-1);
[nr,nphi]=deal(length(r_rg),length(phi_rg));
figure(1);clf;clr=lines();
tlim=5;
[x0,tspan,xuplim]=deal(-sqrt(a)-cm,linspace(-tlim,tlim,2000),6);
[t,x]=deal(cell(1,nphi));
for r=r_rg
    for i=1:nphi
        f=@(t,x)(x+L(t,r,phi_rg(i))).^2-a;
        [t{i},x{i}]=ode45(f,tspan,x0,odeset('RelTol',1e-8,'Events',@(t,x)escape(x,xuplim)));
    end
    tx=[t;x];
    clf;    hold on;
    px=plot(tx{:},'color',clr(1,:),'Linewidth',2,'DisplayName','$x(t)$ ($\dot x=(x+L(t))^2-a$)');
    pe=plot(0,0,'w.','DisplayName','');
    for i=2:nphi
        plot(tspan,L(tspan,r,phi_rg(i)),'Color',(clr(2,:)+1)/2,'LineWidth',1);
    end
    pl=plot(tspan,L(tspan,r,phi_rg(1)),'Color',clr(2,:),'LineWidth',2,...
        'DisplayName','$$L(t)=\frac{4+\frac{1}{2}\sin(\phi-2\pi t)}{1+\exp(-2rt)}$$');
    plot(-tlim*[1,1],sqrt(a)*[-1,1],'o','MarkerFaceColor','k','MarkerSize',10);
    plot(tlim*[1,1],sqrt(a)*[-1,1]-cp,'ko','MarkerFaceColor',0.5*[1,1,1],'MarkerSize',10);
    xlim(tlim*[-1,1]);ylim([-6,6]);
    title(sprintf('a=%g, r=%g',a,r));
    legend([px(1),pe,pl],'Location','northwest','Interpreter','latex','FontSize',20);
    set(gca,'FontSize',18);
    drawnow
    pause
end
%%
function y=Lfun(t,cm,beta,omega,cp,r,phi)
ramp=1./(1+exp(-2*r*t));
y=cm+ramp*(cp-cm)+beta*ramp.*sin(phi-omega*t);
end
%%
function [value,isterminal,direction] = escape(x,xlim)
value=x-xlim;
isterminal=true;
direction=1;
end

