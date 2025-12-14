function [varfields,xtype,kind]=dde_hopf_vars()
varfields=struct('x',1,'v',2,'omega',1,'parameter',1);
kind='stst';
xtype={'x','v'};
end
