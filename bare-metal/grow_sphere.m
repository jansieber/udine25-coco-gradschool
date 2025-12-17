%% continue sphere in 2d
clear
format compact
startup_coco(fullfile(pwd(),'..','coco_2025January28'));
addpath(fullfile(pwd(),'tools'))
%% Cusp normal form
f=@(x,y,z)x.^2+y.^2+z.^2-1;
%% Zero problem Phi1 and initial guess
Phi1=@(u)f(u(1),u(2),u(3));
u0=[1;0;0];
%% Create empty coco problem  and add Phi
prob=coco_prob();
prob=coco_add_func(prob,'sphere',f2coco(Phi1),[],'zero','u0',u0);
%% Add monitor function Psi1 and continuation parameters a,b
% short cut would be 
% prob=coco_add_pars(prob,'',[2,3],{'a','b'});
Psi1=@(u)u;
prob=coco_add_func(prob,'spherepars',f2coco(Psi1),[],'inactive',{'x','y','z'},'uidx',[1,2,3]);
prob=coco_set(prob,'cont','atlas','kd','PtMX',1000,'NAdapt',1,'R',0.1,'R_max',10,'R_min',0.01);
bd=coco(prob,'sphere','',2,{'x','y','z'},{[-2,2],[-2,2],[-2,2]});
%%
atlas=coco_bd_read('sphere','atlas');
figure(1);clf;
plot_atlas_kd(atlas.charts,1,2,3,'polyhedra');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal
[info,atlas]=info_from_run('sphere');
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
axis equal