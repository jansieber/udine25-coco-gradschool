%% Summary
clear
format compact
startup_coco(fullfile(pwd(),'..','coco_2025January28'));
addpath(fullfile(pwd(),'tools'))
%%
u0=[1; 1.1];
prob=coco_prob();
prob=coco_add_func(prob,'Phi1',@circle_J,[],'zero','u0',u0,'f+df');
uidx_Phi1=coco_get_func_data(prob,'Phi1','uidx');
prob=coco_add_func(prob,'Psi1',@Psi1,[],'regular',{'x','y'},'uidx',uidx_Phi1);
prob=coco_set(prob,'cont','PtMX',100);
bd=coco(prob,'circle_demo',[],1,{'x','y'});
figure(1);clf;
coco_plot_bd('circle_demo','x','y');grid on;axis equal
%%  functions for right-hand side for circle with Jacobian provided
function [data,y,J]=circle_J(~,data,u)
y=u(1)^2+(u(2)-1)^2-1;
J=2*[u(1),u(2)-1];
end
%% regular monitoring function Psi1 (no Jacobian needed)
function [data,y]=Psi1(~,data,u)
y=u;
end
