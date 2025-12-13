%% plot composite solution
function plotsol(ax,run,lab,seglist,iv,ip,include_L,include_driver)
if nargin<8
    include_driver=true;
end
if nargin<7
    include_L=false;
end
ramp=@(p,t)1./(1+exp(-2*p(ip.r)*t));
Lfun=@(p,t)p(ip.cm)+ramp(p,t)*(p(ip.cp)-p(ip.cm))+p(ip.beta)*ramp(p,t).*sin(p(ip.phi)-p(ip.omega)*t);
vnames=fieldnames(iv);
segnames={'\gamma','+','-'};
if ~iscell(run)
    run=repmat({run},1,length(seglist));
    lab=repmat({lab},1,length(seglist));
end
for i=length(seglist):-1:1
    sol(i)=bvp_read_solution(seglist{i},run{i},lab{i});
end
t_offset=[sol{2}.T,0,-sol{3}.T];
hw={'off','off','on'};
tbp=@(i)sol{i}.tbp+t_offset(i);
Lt=@(i)Lfun(sol{3}.p,tbp(i));
lw2={'linewidth',2};
lwms=[lw2(:)',{'MarkerSize'},{4}];
lstyle={':','-','-'};
ish=ishold(ax);
hold(ax,'on');
for i=1:length(seglist)
    for k=1:length(vnames)
        if ~include_driver && ~strcmp(vnames{k},'x')
            continue
        end
        plot(ax,sol{i}.tbp+t_offset(i),sol{i}.xbp(:,k),lstyle{i},...
            'DisplayName',['$',vnames{k},'_',segnames{i},'$'],lwms{:});
    end
    if ~include_L
        continue
    end
    plot(ax,tbp(i),-Lt(i),'-','color',0.5*[1,1,1],'HandleVisibility',hw{i},...
        'DisplayName','$-L(t)$',lwms{:});
end
xline(ax,0,'k-','DisplayName','Lin section',lw2{:});
xline(ax,t_offset(1),'k:','DisplayName','$T_+$',lw2{:});
xline(ax,t_offset(3),'k--','DisplayName','$T_-$',lw2{:});
legend(ax,'Interpreter','LaTeX','Location','eastoutside');
title(ax,sprintf('$r=%6.4g$, $a=%6.4g$, $c_+=%4.2g$, $\\beta=%4.2g$',...
    sol{1}.p(ip.r),sol{1}.p(ip.a),sol{1}.p(ip.cp),sol{1}.p(ip.beta)),'Interpreter','latex')
set(ax,'FontSize',18,'box','on','LineWidth',2);
xlabel(ax,'time $t$','Interpreter','latex');
ax.XLim(1)=t_offset(3)*1.05;
ax.YLimMode='auto';
axis(ax,'tight');
%ylim(yrg);
if ~ish
    hold(ax,'off');
end
end