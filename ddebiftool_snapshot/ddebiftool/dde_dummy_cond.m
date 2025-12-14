%% dummy sys_cond condition
function [resi,condi]=dde_dummy_cond(point,pref) %#ok<INUSD>
resi=zeros(0,1);
condi=repmat(point,0,1);
end
