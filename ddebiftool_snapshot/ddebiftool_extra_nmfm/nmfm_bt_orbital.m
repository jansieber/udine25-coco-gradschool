function bt = nmfm_bt_orbital(funcs, bt, varargin)
%% Bogdanov-Takens normal form coefficients
%
% This function assumes that the nullvectors of bt have been computed
% already.
%
% $Id: nmfm_bt.m 309 2018-10-28 19:02:42Z jansieber $
%%
default={'debug',false,'free_pars',[], 'generic_unfolding', true};
options=dde_set_options(default,varargin,'pass_on');
if ~strcmp(bt.kind,'BT')
   display(bt.kind);
   error('NMFM_BT: did not receive a bt point as argument.');
end
Deltas=arrayfun(@(deg){ch_matrix(funcs,point.x,point.parameter,lambda,...
    'deri',deg)},1:5);
[D,dD,ddD,dddD,ddddD]=deal(Deltas{:});

%% normalize the Jordan chains vectors
% q0'q0=1 and q0'*q1=0 from defining system (that is, q1 could be 0)
q0=bt.nvec.q0; 
q1=bt.nvec.q1; 
p0=bt.nvec.p0; 
p1=bt.nvec.p1;

alpha = p0'*dD*q0 +  1/2 * p1'*ddD*q0;
beta= -( p0'*dD*q1 + 1/2*p0'*ddD*q0 + 1/2*p1'*ddD*q1 + 1/6*p1'*dddD*q0 )/alpha;
p0=p0+beta*p1;
p0=p0/alpha;
p1=p1/alpha;
%% Check conditions and normalizations
% note that q1 or p0 may be 0.
if options.debug
    small=@(x)norm(x)<1e-12;
    assert(small(D*q0));          %  A phi0=0
    assert(small(dD*q0+D*q1));    %  A phi1+phi0=0
    assert(small(p1*D));          % phisun1 A = 0
    assert(small(p0*D+p1*dD));    % phisun1 + phisun0 A = 0
    assert(small(q0'*q0-1));      % |phi0(0)|=1 (this is arbitrary)
    assert(small(q0'*q1));        % phi1(0)^T phi0(0)=0 (this is arbitrary)
    assert(small(p0*dD*q0 + 1/2*p1*ddD*q0 - 1)); % <phisun0,phi0> = 1
    assert(small(p1*dD*q1 + 1/2*p1*ddD*q0 - 1)); % <phisun1,phi1> = 1
    assert(small(p0*dD*q1 + 1/2*p0*ddD*q0 + 1/2*p1*ddD*q1 + 1/6*p1*dddD*q0)); % <phisun0,phi1> = 0
end
p1=p1'; p0=p0';

%% abbreviate derivative
F=nmfm_deriv_define(funcs,bt,'free_pars',options.free_pars,varargin{:});
J1 = F.J1;
B  = F.B;
A1 = F.A1;
J2 = F.J2;
C  = F.C;
B1 = F.B1;
A2 = F.A2;
J3 = F.J3;
par0=bt.parameter(:)*0;
dev0=@(v)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))]);
devl=@(v,lambda)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))],'lambda',lambda); %#ok<NASGU>
devlt=@(v)nmfm_dev_fun([v;par0(:,ones(1,size(v,2)))],'lambda',zeros(size(v,2),1),'t',0:size(v,2)-1);
ax=@(a,x) nmfm_dev_ax(a,x);
psi0 = @(h) nmfm_psi0(h.v, p0, p1, Deltas);
psi1 = @(h) nmfm_psi1(h.v, p0, p1, Deltas);

%% define bordered inverse for the characteristic matrix
Bord = [D p1';q0' 0];
take=@(x) x(1:end-1);
Dinv=@(x) take(Bord\[x;0]);

%% auxiliary fucntion to create solution
AINV = @(kappa,w) nmfm_bt_create_solution(kappa,w.v(1:length(q0),:),Deltas,Dinv,devlt);

%% calculate normal form coefficients
phi0 = dev0(q0);
phi1 = devlt([q1,q0]);

a = 1/2*p1*B(phi0,phi0);
b = p0*B(phi0,phi0) + p1*B(phi0,phi1);

h2000hat = AINV(B(phi0,phi0), ax([2*a],[phi1]));
h1100hat = AINV(B(phi0,phi1), ax([b,1],[phi1,h2000hat]));

theta1000 = -1/(12*a)*p1*(3*B(h2000hat,phi0)+C(phi0,phi0,phi0)) + 1/2*psi1(h1100hat);

gamma1 = p0*B(phi0,phi1)-psi0(h2000hat)+1/2*p1*B(phi1,phi1)+theta1000;
gamma2 = 1/(6*a)*(p1*(2*B(phi0,h1100hat)+B(h2000hat,phi1)+C(phi0,phi0,phi1)) + ...
                  2*a*p0*B(phi1,phi1)-b*p1*B(phi1,phi1) + ...
                  p0*(3*B(h2000hat,phi0)+C(phi0,phi0,phi0))+3*gamma1*b - ...
                  10*a*psi0(h1100hat));
h2000 = ax([1;gamma1],[h2000hat,phi0]);
h1100 = ax([1,gamma1-theta1000,gamma2],[h1100hat,phi1,phi0]);
h0200 = AINV(B(phi1,phi1),ax(2,h1100));
h3000 = AINV(3*B(h2000,phi0) + C(phi0,phi0,phi0), ax([6*a,-6*a*theta1000],[h1100,phi1]));
h2100 = AINV(2*B(h1100,phi0) + B(h2000,phi1) + C(phi0,phi0,phi1), ...
            ax([2*a,2*b,1,-2*theta1000*b,2*theta1000^2,-2*theta1000], [h0200,h1100,h3000,phi1,phi0,h2000]));


if options.generic_unfolding
    %% parameter-dependent linear terms H0010,K10,H0001,K01
    nu = (p1*J1)';
    K10hat = nu/(nu'*nu);
    h0010hat = AINV(J1*K10hat,phi1);
    K01hat = [[0,-1];[1,0]]*nu;
    h0001hat = dev0(Dinv(J1*K01hat));

    gamma3 = -(p1*(B(h0001hat,phi0)+A1(phi0,K01hat)))/(2*a);
    delta1 = 1/(p1*(B(h0001hat,phi1)+A1(phi1,K01hat)) + p0*(B(h0001hat,phi0)+A1(phi0,K01hat))+gamma3*b);
    gamma4 = (psi1(h1100)-theta1000-p1*(B(h0010hat,phi0)+A1(phi0,K10hat)))/(2*a);
    delta2 = -p1*(B(h0010hat,phi1)+A1(phi1,K10hat))-gamma4*b+psi1(h0200)-p0*(B(h0010hat,phi0)+A1(phi0,K10hat)) + psi0(h1100);
    K01 = delta1*K01hat;
    h0001 = ax(delta1*[1,gamma3],[h0001hat,phi0]);
    K10 = K10hat+delta2*K01;
    h0010 = ax([1,delta2,gamma4],[h0010hat,h0001,phi0]);

    %% h1010, h0100
    h1010 = AINV(B(h0010,phi0)+A1(phi0,K10), ax([1,-theta1000],[h1100,phi1]));
    h0110 = AINV(B(h0010,phi1)+A1(phi1,K10), ax([1,1],[h0200,h1010]));


    %% h1001, h0101, h2001, h1101
    h1001hat = dev0(Dinv(B(h0001,phi0)+A1(phi0,K01)));
    h0101hat = AINV(B(h0001,phi1)+A1(phi1,K01),ax([1,1],[h1001hat,phi1]));

    %% H1001, H0101, H2001, H1101
    zeta1 = psi1(ax([2*a],[h0101hat])) -p1*(A1(h2000,K01)+B(h0001,h2000)+2*B(h1001hat,phi0) ...
                    +B1(phi0,phi0,K01)+C(h0001,phi0,phi0));
    zeta2 = psi1(ax([b,1,-theta1000],[h0101hat,h1100,h1001hat])) - theta1000 -p1*(A1(h1100,K01) ...
                    +B(h0001,h1100)+B(h0101hat,phi0) ...
                    +B(h1001hat,phi1)+B1(phi0,phi1,K01)+C(h0001,phi0,phi1)) ...
                    +psi0(ax(2*a,h0101hat)) - p0*(A1(h2000,K01)+B(h0001,h2000) ...
                    +2*B(h1001hat,phi0)+B1(phi0,phi0,K01)+C(h0001,phi0,phi0));
    gamma5 = 2*zeta2/b-zeta1/(2*a);
    theta0001 = zeta1/(2*a)-zeta2/b;
    h1001 = ax([1,gamma5],[h1001hat,phi0]);
    h0101 = ax([1,gamma5,-theta0001],[h0101hat,phi1,phi1]);
    h2001 = AINV(+A1(h2000,K01)+B(h0001,h2000)+2*B(h1001,phi0) ...
                    +B1(phi0,phi0,K01)+C(h0001,phi0,phi0),ax(2*a*[1,-theta0001],[h0101,phi1]));
    h1101 = AINV(A1(h1100,K01) + B(h0001,h1100) +B(h0101,phi0)+ B(h1001,phi1) +B1(phi0,phi1,K01)+C(h0001,phi0,phi1), ...
                    ax([b,1,1,-theta1000*[1,1,-theta0001],-theta0001*[1,b,-theta1000]], ...
                    [h0101,h1100,h2001,h1001,phi1,phi0,h2000,phi1,phi0]));

    %% K11, h0011
    K11 = (psi1(h0101) - theta0001 - p1*(A1(h0001,K10)+A1(h0010,K01)+B(h0001,h0010)+J2(K01,K10)))*K10;
    h0011 = AINV(J1*K11+A1(h0001,K10)+A1(h0010,K01)+B(h0001,h0010)+J2(K01,K10),ax([1,-theta0001],[h0101,phi1]));

    %% K02, h0002, h1002, h0102
    K02hat = -(p1*(2*A1(h0001,K01) + B(h0001,h0001) + J2(K01,K01)))*K10;

    h0002hat = dev0(Dinv(J1*K02hat + 2*A1(h0001,K01) + B(h0001,h0001) + J2(K01,K01)));
    gamma6 = -(1/(2*a))*p1*(2*A1(h1001,K01)+A1(phi0,K02hat)+A2(phi0,K01,K01) ...
                    +B(phi0,h0002hat)+2*B(h0001,h1001)+2*B1(phi0,h0001,K01) ...
                    +C(phi0,h0001,h0001));
    delta3 = -psi1(ax(2*[theta0001,-1], [h1001,h0101])) - 2*theta0001  ...
                    -p1*(2*A1(h0101,K01)+A1(phi1,K02hat)+A2(phi1,K01,K01) ...
                    +B(phi1,h0002hat)+2*B(h0001,h0101)+2*B1(phi1,h0001,K01) ...
                    +C(phi1,h0001,h0001)) ...
                    -p0*(2*A1(h1001,K01)+A1(phi0,K02hat)+A2(phi0,K01,K01) ...
                    +B(phi0,h0002hat)+2*B(h0001,h1001)+2*B1(phi0,h0001,K01) ...
                    +C(phi0,h0001,h0001))-gamma6*b;
    K02 = K02hat + delta3*K01;
    h0002 = ax([1,delta3,gamma6],[h0002hat,h0001,phi0]);
    h1002 = dev0(Dinv(2*A1(h1001,K01)+A1(phi0,K02)+A2(phi0,K01,K01)+B(phi0,h0002) ...
                    +2*B(h0001,h1001)+2*B1(phi0,h0001,K01)+C(phi0,h0001,h0001)));
    h0102 = AINV(2*A1(h0101,K01)+A1(phi1,K02)+A2(phi1,K01,K01)+B(phi1,h0002) ...
                    +2*B(h0001,h0101)+2*B1(phi1,h0001,K01)+C(phi1,h0001,h0001), ...
                    ax([2,1,-2*theta0001*[1,1,-theta0001]],[h0101,h1002,h1001,phi1,phi0]));

    %% K03, h0003
    K03 = -p1*(A1(h0001,K02)+A1(h0002,K01)+2*(A1(h0001,K02) ...
                    +A1(h0002,K01))+3*B(h0001,h0002)+3*J2(K01,K02) ...
                    +3*A2(h0001,K01,K01)+3*B1(h0001,h0001,K01)+C(h0001,h0001,h0001) ...
                    +J3(K01,K01,K01))*K10;
    h0003 = dev0(Dinv(J1*K03+A1(h0001,K02)+A1(h0002,K01)+2*(A1(h0001,K02) ...
                    +A1(h0002,K01))+3*B(h0001,h0002)+3*J2(K01,K02) ...
                    +3*A2(h0001,K01,K01)+3*B1(h0001,h0001,K01)+C(h0001,h0001,h0001) ...
                    +J3(K01,K01,K01)));

    bt.nmfm.a2=a;
    bt.nmfm.b2=b;
    bt.nvec.p1=p1';
    bt.nvec.p0=p0';
    % bt.nmfm.a3=1/6*p1*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0,PHI_0);
    % bt.nmfm.b3=1/2*p1*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0,PHI_1) ...
    %     + 1/2*p0*sys_mfderi(xx,bt.parameter,PHI_0,PHI_0,PHI_0);

    bt.nmfm.a=a;
    bt.nmfm.b=b;
    bt.nmfm.theta1000 = theta1000;
    bt.nmfm.theta0001 = theta0001;
    bt.nmfm.K10 = K10;
    bt.nmfm.K01 = K01;
    bt.nmfm.K02 = K02;
    bt.nmfm.K11 = K11;
    bt.nmfm.K03 = K03;
    bt.nmfm.phi0 = phi0;
    bt.nmfm.phi1 = phi1;
    bt.nmfm.h0010 = h0010;
    bt.nmfm.h0001 = h0001;
    bt.nmfm.h2000 = h2000;
    bt.nmfm.h1100 = h1100;
    bt.nmfm.h0200 = h0200;
    bt.nmfm.h1010 = h1010;
    bt.nmfm.h1001 = h1001;
    bt.nmfm.h0110 = h0110;
    bt.nmfm.h0101 = h0101;
    bt.nmfm.h0002 = h0002;
    bt.nmfm.h0011 = h0011;
    bt.nmfm.h3000 = h3000;
    bt.nmfm.h2100 = h2100;
    bt.nmfm.h1101 = h1101;
    bt.nmfm.h2001 = h2001;
    bt.nmfm.h0003 = h0003;
    bt.nmfm.h1002 = h1002;
    bt.nmfm.h0102 = h0102;
else
    e1 = [1;0];
    e2 = [0;1];
    M0 = [        p1*A1(phi0,e1)                     p1*A1(phi0,e2) 
     (p0*A1(phi0,e1) + p1*A1(phi1,e1))     (p0*A1(phi0,e2) + p1*A1(phi1,e2)) ];
    Id = [1 0; 0 1];
    delta = M0\Id;

    K10 = delta(1)*e1 + delta(2)*e2;
    K01 = delta(3)*e1 + delta(4)*e2;

    h1010hat = AINV(A1(phi0,K10), phi1);
    h1001hat = dev0(Dinv(A1(phi0,K01)));

    h0110hat = AINV(A1(phi1,K10),h1010hat);
    h0101hat = AINV(A1(phi1,K01), ax([1, 1], [h1001hat, phi1]));

    zeta1 = -p1*(A1(h2000,K10) + 2*B(h1010hat,phi0) + B1(phi0,phi0,K10)) + 2*psi1(ax([a,1],[h0110hat, h1100])) - 2*theta1000;
    zeta2 = -p0*(A1(h2000,K10) + 2*B(h1010hat,phi0) + B1(phi0,phi0,K10)) + 2*psi0(ax([a,1],[h0110hat, h1100])) - ...
            p1*(A1(h1100,K10) + B(h0110hat,phi0) + B(h1010hat,phi1) + B1(phi0,phi1,K10)) +  ...
            psi1(ax([b,1,-theta1000],[h0110hat, h0200, h1010hat]));
    gamma3 = 2*zeta2/b-zeta1/(2*a);
    theta0010 = zeta1/(2*a)-zeta2/b;

    % gamma3 = zeta1/(2*a)
    % theta0010 = 0

    zeta3 = -p1*(A1(h2000,K01) + 2*B(h1001hat,phi0) + B1(phi0,phi0,K01)) + 2*a*psi1(h0101hat);
    zeta4 =  p0*(A1(h2000,K01) + 2*B(h1001hat,phi0) + B1(phi0,phi0,K01)) + 2*a*psi0(h0101hat) - ... 
            p1*(A1(h1100,K01) + B(h0101hat,phi0) + B(h1001hat,phi1) + B1(phi0,phi1,K01)) +  ...
            psi1(ax([b,1,-theta1000], [h0101hat, h1100, h1001hat])) + theta1000;
    gamma4 = 2*zeta4/b-zeta3/(2*a);
    theta0001 = zeta3/(2*a)-zeta4/b;

    % gamma4 = zeta3/(2*a)
    % theta0001 = 0

    gamma5 = -p1*(2*A1(h1001hat,K01) + A2(phi0,K01,K01));
    gamma6 = -p0*(2*A1(h1001hat,K01) + A2(phi0,K01,K01)) -p1*(2*A1(h0101hat,K01) + A2(phi1,K01,K01)) + 2*psi1(h0101hat) - 2*theta0001;

    gamma7 = -p1*(A1(h1001hat,K10) + A1(h1010hat,K01) + A2(phi0,K01,K10)) + psi1(h0101hat) + - 2*theta0001;
    gamma8 = -p0*(A1(h1001hat,K10) + A1(h1010hat,K01) + A2(phi0,K01,K10)) - ...
          p1*(A1(h0101hat,K10) + A1(h0110hat,K01) + A2(phi1,K01,K10)) - theta0010 + psi0(h0101hat) + psi1(h0110hat);

    gamma9  = -p1*(2*A1(h1010hat,K10) + A2(phi0,K10,K10)) - 4*theta0010 + 2*psi1(h0110hat);
    gamma10 = -p0*(2*A1(h1010hat,K10) + A2(phi0,K10,K10)) - p1*(2*A1(h0110hat,K10) + A2(phi1,K10,K10)) + 2*psi0(h0110hat);

    h1010 = ax([1,gamma3], [h1010hat, phi0]);
    h1001 = ax([1,gamma4], [h1001hat, phi0]);
    h0110 = ax([1,gamma3,-theta0010], [h0110hat, phi1, phi1]);
    h0101 = ax([1,gamma4,-theta0001], [h0101hat, phi1, phi1]);

    K02 = gamma5*K10 +  gamma6*K01;
    K11 = gamma7*K10 +  gamma8*K01;
    K20 = gamma9*K10 + gamma10*K01;

    h2010 = AINV(A1(h2000,K10) + 2*B(h1010,phi0) + B1(phi0,phi0,K10), ax(2*[a, 1, -theta0010*a, -theta1000], [h0110, h1100, phi1, phi1]));
    h1110 = AINV(A1(h1100,K10) + B(h0110,phi0) + B(h1010,phi1) + B1(phi0,phi1,K10), ax([b,1,1,-theta1000*[1,-theta0010],-theta0010*[b,-theta1000,1]], [h0110, h0200, h2010, h1010, phi0, phi1, phi0, h2000]));

    h2001 = AINV(A1(h2000,K01) + 2*B(h1001,phi0) + B1(phi0,phi0,K01), ax(2*a*[1,-theta0001], [h0101, phi1]));
    h1101 = AINV(A1(h1100,K01) + B(h0101,phi0) + B(h1001,phi1) + B1(phi0,phi1,K01), ax([b,1,1,-theta1000*[1,1,-theta0001],-theta0001*[b,-theta1000,1]], [h0101, h1100, h2001, h1001, phi1, phi0, phi1, phi0, h2000]));

    h1002 = dev0(Dinv( 2*A1(h1001,K01) + A1(phi0,K02) + A2(phi0,K01,K01)));
    h0102 = AINV(2*A1(h0101,K01) + A1(phi1,K02) + A2(phi1,K01,K01), ax([2,1,-2*theta0001*[1,1,-theta0001]], [h0101, h1002, h1001, phi1, phi0]));

    h1020 = AINV(2*A1(h1010,K10) + A1(phi0,K20) + A2(phi0,K10,K10), ax(2*[1, -theta0010], [h0110, phi1]));
    h0120 = AINV(2*A1(h0110,K10) + A1(phi1,K20) + A2(phi1,K10,K10), ax([1,-2*theta0010*[1,-theta0010]], [h1020, h1010, phi0]));

    h1011 = AINV(A1(h1001,K10) + A1(h1010,K01) + A1(phi0,K11) + A2(phi0,K01,K10), ax([1,-theta0001], [h0101, phi1]));
    h0111 = AINV(A1(h0101,K10) + A1(h0110,K01) + A1(phi1,K11) + A2(phi1,K01,K10), ax([1,1,-theta0010*[1,1,-theta0001],-theta0001*[1,-theta0010]], [h0110, h1011, h1001, phi1, phi0, h1010, phi0]));

    bt.nmfm.a = a; 
    bt.nmfm.b = b;
    bt.nmfm.theta1000 = theta1000;
    bt.nmfm.theta0010 = theta0010;
    bt.nmfm.theta0001 = theta0001;
    bt.nmfm.phi0    = phi0;
    bt.nmfm.phi1    = phi1;
    bt.nmfm.h2000 = h2000;
    bt.nmfm.h1100 = h1100;
    bt.nmfm.h0200 = h0200;
    bt.nmfm.h3000 = h3000;
    bt.nmfm.h2100 = h2100;
    bt.nmfm.K10   = K10;
    bt.nmfm.K01   = K01;
    bt.nmfm.K02   = K02;
    bt.nmfm.K11   = K11;
    bt.nmfm.K20   = K20;
    bt.nmfm.h1010 = h1010;
    bt.nmfm.h1001 = h1001;
    bt.nmfm.h0110 = h0110;
    bt.nmfm.h0101 = h0101;
    bt.nmfm.h2010 = h2010;
    bt.nmfm.h1110 = h1110;
    bt.nmfm.h2001 = h2001;
    bt.nmfm.h1101 = h1101;
    bt.nmfm.h1002 = h1002;
    bt.nmfm.h0102 = h0102;
    bt.nmfm.h1020 = h1020;
    bt.nmfm.h0120 = h0120;
    bt.nmfm.h1011 = h1011;
    bt.nmfm.h0111 = h0111;
end

n = length(bt.x);
take = @(v) v(1:n);
% bt = BT_nmfm_orbital(odefile, bt, ap, options);
phi0    = take(nmfm_dev_call(bt.nmfm.phi0,0));
phi1    = take(nmfm_dev_call(bt.nmfm.phi1,0));
h2000 = take(nmfm_dev_call(bt.nmfm.h2000,0));
h1100 = take(nmfm_dev_call(bt.nmfm.h1100,0));
h0200 = take(nmfm_dev_call(bt.nmfm.h0200,0));
h3000 = take(nmfm_dev_call(bt.nmfm.h3000,0));
h2100 = take(nmfm_dev_call(bt.nmfm.h2100,0));
K10   = bt.nmfm.K10;
K01   = bt.nmfm.K01;
K02   = bt.nmfm.K02;
K11   = bt.nmfm.K11;
h1010 = take(nmfm_dev_call(bt.nmfm.h1010,0));
h1001 = take(nmfm_dev_call(bt.nmfm.h1001,0));
h0110 = take(nmfm_dev_call(bt.nmfm.h0110,0));
h0101 = take(nmfm_dev_call(bt.nmfm.h0101,0));
h1101 = take(nmfm_dev_call(bt.nmfm.h1101,0));
h2001 = take(nmfm_dev_call(bt.nmfm.h2001,0));
h1002 = take(nmfm_dev_call(bt.nmfm.h1002,0));
h0102 = take(nmfm_dev_call(bt.nmfm.h0102,0));
if options.generic_unfolding
    K03   = bt.nmfm.K03;
    h0010 = take(nmfm_dev_call(bt.nmfm.h0010,0));
    h0001 = take(nmfm_dev_call(bt.nmfm.h0001,0));
    h0002 = take(nmfm_dev_call(bt.nmfm.h0002,0));
    h0011 = take(nmfm_dev_call(bt.nmfm.h0011,0));
    h0003 = take(nmfm_dev_call(bt.nmfm.h0003,0));
    bt.nmfm.K = @(beta1,beta2) K10*beta1 + K01*beta2 + 1/2*K02*beta2.^2 + K11*beta1.*beta2 + 1/6*K03*beta2^3;
    bt.nmfm.H = @(w0,w1,beta1,beta2) phi0.*w0 + phi1.*w1 + h0010*beta1 + h0001*beta2 + ...
                    1/2*h2000*w0.^2 + h1100*w0.*w1 + 1/2*h0200*w1.^2 +...
                    h1010*w0.*beta1 + h1001*w0.*beta2 + h0110*w1.*beta1 + ...
                    h0101*w1.*beta2 + 1/2*h0002*beta2.^2 + h0011*beta1.*beta2 + ...
                    1/6*h3000*w0.^3 + 1/2*h2100*w0.^2.*w1 + h1101*w0.*w1.*beta2 + ...
                    1/2*h2001*w0.^2*beta2 + 1/6*h0003*beta2.^3 + ...
                    1/2*h1002*w0*beta2.^2 + 1/2*h0102*w1*beta2.^2;
else
    K20 = bt.nmfm.K20;
    h2010 = take(nmfm_dev_call(bt.nmfm.h2010,0));
    h1110 = take(nmfm_dev_call(bt.nmfm.h1110,0));
    h1020 = take(nmfm_dev_call(bt.nmfm.h1020,0));
    h0120 = take(nmfm_dev_call(bt.nmfm.h0120,0));
    h1011 = take(nmfm_dev_call(bt.nmfm.h1011,0));
    h0111 = take(nmfm_dev_call(bt.nmfm.h0111,0));
    bt.nmfm.K = @(beta1,beta2) K10*beta1 + K01*beta2 + 1/2*K20*beta1^2 + K11*beta1*beta2 + 1/2*K02*beta2^2;
    bt.nmfm.H = @(w0,w1,beta1,beta2) phi0.*w0 + phi1.*w1 + 1/2.*h2000.*w0.^2 + h1100.*w0.*w1 + 1/2.*h0200.*w1.^2 + h1010.*w0.*beta1 + ...
                    h1001.*w0.*beta2 + h0110.*w1.*beta1 + h0101.*w1.*beta2 + 1/2.*h0102.*w1.*beta2.^2 + ...
                    h0111.*w1.*beta1.*beta2 + 1/2.*h0120.*w1.*beta1.^2 + 1/2.*h1002.*w0.*beta2.^2 + h1011.*w0.*beta1.*beta2 + ...
                    1/2.*h1020.*w0.*beta1.^2 + h1101.*w0.*w1.*beta2 + h1110.*w0.*w1.*beta1 +  ...
                    1/2.*h2001.*w0.^2.*beta2 + 1/2.*h2010.*w0.^2.*beta1 + 1/2.*h2100.*w0.^2.*w1 + 1/6.*h3000.*w0.^3;
end

end
