%% Illustration of adding continuation parameters for kd continuation
% execute coco's startup for setting path
clear
format compact
%% Cusp normal form
f=@(x,a,b)a+b*x-x^3;
%% Zero problem Phi1 and initial guess
Phi1=@(u)f(u(1),u(2),u(3));
u0=[1;0;1];
%% Create empty coco problem  and add Phi
prob=coco_prob();
prob=coco_add_func(prob,'cusp',f2coco(Phi1),[],'zero','u0',u0);
%% Add monitor function Psi1 and continuation parameters a,b
% short cut would be 
% prob=coco_add_pars(prob,'',[2,3],{'a','b'});
Psi1=@(u)u;
prob=coco_add_func(prob,'cusppars',f2coco(Psi1),[],'inactive',{'x','a','b'},'uidx',[1,2,3]);
prob=coco_set(prob,'cont','atlas','kd','PtMX',1000,'NAdapt',1,'R',0.1,'R_max',10,'R_min',0.01);
bd=coco(prob,'arun','',2,{'a','x','b'},{[-1,2],[-3,3],[-5,5]});
%%
atlas=coco_bd_read('arun','atlas');
figure(1);clf;
plot_atlas_kd(atlas.charts,1,2,3,'polyhedra');
grid on;
xlabel('x');
ylabel('a');
zlabel('b');
set(gca,'FontSize',18);
figure(2);
[info,atlas]=info_from_run('arun');
xp=info.xp;
adj=adjacency_from_atlas(atlas);
tri=triangulate_adjacency(adj);
figure(3);clf;
trisurf(tri,xp(:,1),xp(:,2),xp(:,3),'edgecolor','k');
grid on;
xlabel(info.pnames{1});
ylabel(info.pnames{2});
zlabel(info.pnames{3});
set(gca,'FontSize',18);
%% Add further monitor function Psi2 and continuation parameters a,b,c 
% to stay close to saddle-node
% prob=coco_add_pars(prob,'',[2,3],{'a','b'});
probs=coco_prob();
probs=coco_add_func(probs,'cusp',f2coco(Phi1),[],'zero','u0',[0;0;0]);
Psi2=@(u)u(2)-3*u(1)^2;
probs=coco_add_pars(probs,'cusppars',[1,2,3],{'x','a','b'});
probs=coco_add_func(probs,'sndiff',f2coco(Psi2),[],'active',{'c'},'uidx',[1,3]);
probs=coco_set(probs,'cont','atlas','kd','PtMX',1000,...
    'NAdapt',1,'R',0.1,'R_max',1,'R_min',0.05);
bds=coco(probs,'snrun','',2,{'a','x','b','c'},{[-1,2],[-3,3],[-2,2],[-0.5,0.5]});
%%
sn_atlas=coco_bd_read('snrun','atlas');
figure(2);clf;
plot_atlas_kd(sn_atlas.charts,3,'polyhedra');
%plot_atlas_kd(sn_atlas.charts,3,'basepoints');
grid on;
xlabel('x');
ylabel('a');
zlabel('b');
%%
[info,sn_atlas]=info_from_run('snrun');
xp=info.xp;
adj=adjacency_from_atlas(sn_atlas);
tri=triangulate_adjacency(adj);
figure(3);clf;
trisurf(tri,xp(:,1),xp(:,2),xp(:,3),'edgecolor','k');
grid on;
xlabel(info.pnames{1});
ylabel(info.pnames{2});
zlabel(info.pnames{3});
