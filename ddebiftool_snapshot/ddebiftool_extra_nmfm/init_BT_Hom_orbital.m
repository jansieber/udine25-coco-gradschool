function [hcli,ups, p] = init_BT_Hom_orbital(funcs,bt,varargin) 

default = {'debug',0,'free_pars',[],'ntst',40,'ncol',4,'order',3, ...
    'TTolerance',zeros(2,1),'amplitude',zeros(2,1),'method','orbital', ...
    'perturbationparameter',[0.1; 0.1] ...
    'amplitudeTToleranceRatio',[1.0e-03; 1.0e-03],'generic_unfolding',true, ...
    'degree',3};
options = dde_set_options(default,varargin,'pass_on');

%% Set coefficients
n = length(bt.x);
take = @(v) v(1:n);
% bt = BT_nmfm_orbital(odefile, bt, ap, options);
a  = bt.nmfm.a;
b  = bt.nmfm.b;
theta1000 = bt.nmfm.theta1000;
theta0001 = bt.nmfm.theta0001;
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
else
    theta0010 = bt.nmfm.theta0010;
    K20 = bt.nmfm.K20;
    h2010 = take(nmfm_dev_call(bt.nmfm.h2010,0));
    h1110 = take(nmfm_dev_call(bt.nmfm.h1110,0));
    h1020 = take(nmfm_dev_call(bt.nmfm.h1020,0));
    h0120 = take(nmfm_dev_call(bt.nmfm.h0120,0));
    h1011 = take(nmfm_dev_call(bt.nmfm.h1011,0));
    h0111 = take(nmfm_dev_call(bt.nmfm.h0111,0));
end

if options.debug
    disp('BT normal form coefficients:')
    fprintf('a = %d,\t b=%d\n', a, b);
end

if options.generic_unfolding
    number_of_orbits = 1;
else
    hcli = cell(2,1);
    number_of_orbits = 2;
    if isscalar(options.perturbationparameter)
        options.perturbationparameter = options.perturbationparameter*[1,1];
    end
end

for i=1:number_of_orbits
    %% Set amplitude and TTolerance
    if ~isempty(options.amplitude(i)) && options.amplitude(i) ~= 0
        amplitude(i) = options.amplitude(i);
        eps = sqrt(amplitude(i)*b^2/(6*abs(a)));
    else
        eps = options.perturbationparameter(i);
        amplitude = eps^2/b^2*6*abs(a);
    end
    if ~isempty(options.TTolerance(i)) && options.TTolerance(i) ~= 0
        TTolerance = options.TTolerance(i);
    else
        TTolerance = amplitude*options.amplitudeTToleranceRatio(i);
    end
    if options.debug
        fprintf('The initial perturbation parameter epsilon:  %d\n', eps);
        fprintf('The initial amplitude: %g\n', amplitude);
    end

    %% The initial approximation of cycle
    % transformation of time \xi to s
    xi1 = @(s) -6.*log(cosh(s))/7;
    xi2 = @(s) (-9/98).*(4.*s+(-5).*tanh(s)+(-8).*log( cosh(s)).*tanh(s));
    xi3 = @(s) (9/2401).*((-91)+(-92).*log(cosh(s))+(91+21.*(3+(-4).*log(cosh(s)) ...
       ).*log(cosh(s))).*sech(s).^2+84.*s.*tanh(s));

    % homoclinic orbit
    u0 = @(xi) 6*tanh(xi).^2 - 4; 
    u2 = @(xi) -18/49*sech(xi).^2; 
    v0 = @(xi) 12.*sech(xi).^2.*tanh(xi);
    v1 = @(xi) -72/7*tanh(xi).*sech(xi).^2.*tanh(xi);
    v2 = @(xi) (90/49 + 162/49*tanh(xi).^2).*sech(xi).^2.*tanh(xi);
    v3 = @(xi) -(3888/2401*tanh(xi) - 216/343*tanh(xi).^3).*sech(xi).^2.*tanh(xi);
    switch options.order
        case 0
            xi  = @(s) s;
            u   = @(xi) u0(xi);
            v   = @(xi) v0(xi);
        case 1
            xi  = @(s) s + xi1(s).*eps;
            u   = @(xi) u0(xi);
            v   = @(xi) v0(xi)+v1(xi).*eps;
        case 2
            xi  = @(s) s + xi1(s).*eps + xi2(s).*eps^2;
            u   = @(xi) u0(xi)+u2(xi).*eps.^2;
            v   = @(xi) v0(xi)+v1(xi).*eps+v2(xi).*eps.^2;
        otherwise
            xi  = @(s) s + xi1(s).*eps + xi2(s).*eps^2 + xi3(s).*eps^3;
            u   = @(xi) u0(xi)+u2(xi).*eps.^2;
            v   = @(xi) v0(xi)+v1(xi).*eps+v2(xi).*eps.^2+v3(xi).*eps.^3;
    end
    if options.generic_unfolding
        w0 = @(tau) a/b^2*u(xi(a/b*eps*tau))*eps^2;
    else
        if i==1
            w0 = @(tau) a/b^2*(u(xi(a/b*eps*tau))-2)*eps^2;
        else
            w0 = @(tau) a/b^2*(u(xi(a/b*eps*tau))+2)*eps^2;
        end
    end
    w1 = @(tau) a^2/b^3*v(xi(a/b*eps*tau))*eps^3;

    %% Approximate parameters
    tau0  = 10/7;
    tau2  = 288/2401;

    if options.generic_unfolding
        beta1 = -4*a^3/b^4*eps^4;
        switch options.order
            case 1
                beta2 = a/b*tau0*eps^2;
                alpha = K10*beta1 + K01*beta2 + 1/2*K02*beta2.^2;
            otherwise
                beta2 = a/b*(tau0 + tau2*eps^2)*eps^2;
                alpha = K10*beta1 + K01*beta2 + 1/2*K02*beta2.^2 + K11*beta1.*beta2 + 1/6*K03*beta2^3;
        end
    else
        if i==1
            beta1 =  4*a^2/b^2*eps^2;
        else
            beta1 = -4*a^2/b^2*eps^2;
        end
        switch options.order
            case 1
                tau = tau0;
                if i==1
                    beta2 = a/b*(tau + 2)*eps^2;
                else
                    beta2 = a/b*(tau - 2)*eps^2;
                end
                alpha = K10.*beta1 + K01.*beta2;
            otherwise
                tau = tau0 + tau2*eps^2;
                if i==1
                    beta2 = a/b*(tau + 2)*eps^2;
                else
                    beta2 = a/b*(tau - 2)*eps^2;
                end
                alpha = K10.*beta1 + K01.*beta2 + 1/2*K20.*beta1.^2 + K11.*beta1.*beta2 + 1/2*K02.*beta2.^2;
        end
    end
    p = bt.parameter(options.free_pars);
    p(options.free_pars) = p(options.free_pars) + alpha';

    %% tau of t
    switch options.order
        case 1
            auxilaryIntegral = @(xi) ...
                2.*xi-6*tanh(xi)+eps.*((-18/7)+(12/7).*log(cosh(xi))+(18/7).*sech(xi).^2);

            tOfTau = @(tau) (1+(10/7).*a.*b.^(-1).*eps.^2.*theta0001).*tau ...
                        + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps));
        case 2
            auxilaryIntegral = @(xi) ... 
                2.*xi+eps.*((-18/7)+(12/7).*log(cosh(xi))+(18/7).*sech(xi).^2)+( ...
                -6).*tanh(xi)+eps.^2.*((36/49).*xi+(-81/49).*tanh(xi)+(45/49).* ...
                sech(xi).^2.*tanh(xi));

            tOfTau = @(tau) (1+a/b*(10/7+288/2401*eps^2)*eps.^2*theta0001)*tau ...
                        + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps));
        otherwise
            auxilaryIntegral = @(xi) ...
                2.*xi+eps.*((-18/7)+(12/7).*log(cosh(xi))+(18/7).*sech(xi).^2)+ ...
                eps.^3.*((-468/2401)+(144/2401).*log(cosh(xi))+(846/2401).*sech( ...
                xi).^2+(-54/343).*sech(xi).^4)+(-6).*tanh(xi)+eps.^2.*((36/49).* ...
                xi+(-81/49).*tanh(xi)+(45/49).*sech(xi).^2.*tanh(xi));

            tOfTau = @(tau) (1+a/b*(10/7+288/2401*eps^2)*eps.^2*theta0001)*tau ...
                        + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps));
    end
    if ~options.generic_unfolding
        switch options.order
            case 1
                if i==1
                    tOfTau = @(tau) tau + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps)) - 2*theta1000*a/b^2*eps^2*tau;
                else
                    tOfTau = @(tau) tau + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps)) + 2*theta1000*a/b^2*eps^2*tau;
                end
            otherwise
                if i==1
                    tOfTau = @(tau) (1 + beta1*theta0010 + beta2*theta0001)*tau + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps)) - 2*theta1000*a/b^2*eps^2*tau;
                else
                    tOfTau = @(tau) (1 + beta1*theta0010 + beta2*theta0001)*tau + theta1000*1/b*eps*auxilaryIntegral(xi(a/b*tau*eps)) + 2*theta1000*a/b^2*eps^2*tau;
                end
            end
    end

    % Estimate half-return time 
    tauofTplus = abs(b/a*1/eps*asech(abs(b)/eps*sqrt(TTolerance/(6*abs(a)))));
    T = fzero(@(tau) tauofTplus - tOfTau(tau), tauofTplus);

    if options.debug
        fprintf('The initial half-return time T: %g\n', T);
    end
    dt = 1/(options.ntst*options.ncol-1);
    t = 0:dt:1;
    uniformmesh = t;
    t = (2*t - 1)*T;
    tauOfT = zeros(size(t));
    for j = 1:length(t)
        tauOfT(j) = fzero(@(tau) t(j) - tOfTau(tau), t(j));
    end
    w0 = w0(tauOfT);
    w1 = w1(tauOfT);

    %% lift the orbit
    if options.generic_unfolding
        switch options.order
            case 1
                H = @(w0,w1,beta1,beta2) phi0*w0 + phi1*w1 + h0010*beta1 + h0001*beta2 + ...
                    1/2*h2000*w0.^2 + h1100*w0.*w1 + ...
                    h1001*w0.*beta2 + h0101*w1.*beta2 + 1/2*h0002*beta2.^2;
            case 2
                H = @(w0,w1,beta1,beta2) phi0*w0 + phi1*w1 + h0010*beta1 + h0001*beta2 + ...
                    1/2*h2000*w0.^2 + h1100*w0.*w1 + 1/2*h0200*w1.^2 +...
                    h1010*w0.*beta1 + h1001*w0.*beta2 +  ...
                    h0101*w1.*beta2 + 1/2*h0002*beta2.^2 + h0011*beta1.*beta2 + ...
                    1/6*h3000*w0.^3 + ...
                    1/2*h2001*w0.^2*beta2 + 1/6*h0003*beta2.^3 + ...
                    1/2*h1002*w0*beta2.^2;
            otherwise
                H = @(w0,w1,beta1,beta2) phi0*w0 + phi1*w1 + h0010*beta1 + h0001*beta2 + ...
                    1/2*h2000*w0.^2 + h1100*w0.*w1 + 1/2*h0200*w1.^2 +...
                    h1010*w0.*beta1 + h1001*w0.*beta2 + h0110*w1.*beta1 + ...
                    h0101*w1.*beta2 + 1/2*h0002*beta2.^2 + h0011*beta1.*beta2 + ...
                    1/6*h3000*w0.^3 + 1/2*h2100*w0.^2.*w1 + h1101*w0.*w1.*beta2 + ...
                    1/2*h2001*w0.^2*beta2 + 1/6*h0003*beta2.^3 + ...
                    1/2*h1002*w0*beta2.^2 + 1/2*h0102*w1*beta2.^2;
        end
    else
        switch options.order
            case 1
                H = @(w0,w1,beta1,beta2) phi0.*w0 + phi1.*w1 + 1/2.*h2000.*w0.^2 + h1100.*w0.*w1 + h1010.*w0.*beta1 + ...
                                h1001.*w0.*beta2 + h0110.*w1.*beta1 + h0101.*w1.*beta2;
            case 2
                H = @(w0,w1,beta1,beta2) phi0.*w0 + phi1.*w1 + 1/2.*h2000.*w0.^2 + h1100.*w0.*w1 + 1/2.*h0200.*w1.^2 + h1010.*w0.*beta1 + ...
                                h1001.*w0.*beta2 + h0110.*w1.*beta1 + h0101.*w1.*beta2 + ...
                                1/2.*h1002.*w0.*beta2.^2 + h1011.*w0.*beta1.*beta2 + ...
                                1/2.*h1020.*w0.*beta1.^2 +  ...
                                1/2.*h2001.*w0.^2.*beta2 + 1/2.*h2010.*w0.^2.*beta1 + 1/6.*h3000.*w0.^3;
            otherwise
                H = @(w0,w1,beta1,beta2) phi0.*w0 + phi1.*w1 + 1/2.*h2000.*w0.^2 + h1100.*w0.*w1 + 1/2.*h0200.*w1.^2 + h1010.*w0.*beta1 + ...
                                h1001.*w0.*beta2 + h0110.*w1.*beta1 + h0101.*w1.*beta2 + 1/2.*h0102.*w1.*beta2.^2 + ...
                                h0111.*w1.*beta1.*beta2 + 1/2.*h0120.*w1.*beta1.^2 + 1/2.*h1002.*w0.*beta2.^2 + h1011.*w0.*beta1.*beta2 + ...
                                1/2.*h1020.*w0.*beta1.^2 + h1101.*w0.*w1.*beta2 + h1110.*w0.*w1.*beta1 +  ...
                                1/2.*h2001.*w0.^2.*beta2 + 1/2.*h2010.*w0.^2.*beta1 + 1/2.*h2100.*w0.^2.*w1 + 1/6.*h3000.*w0.^3;
        end
    end
    ups = bt.x + H(w0,w1,beta1,beta2);

    %% Approximate the sadddle equilibrium
    delta0 = 2;
    delta2 = 0;
    u0inf = delta0;
    u2inf = delta2;
    uinf  = u0inf+u2inf.*eps.^2;
    if options.generic_unfolding
        w0inf = a/b^2*uinf.*eps.^2;
    else
        if i==1
            uinf = uinf - 2;
        else
            uinf = uinf + 2;
        end
        w0inf = a/b^2*uinf.*eps.^2;
    end
    x0 = bt.x + H(w0inf,0,beta1,beta2);

    %% construct homoclinic structure
    hcli{i}.kind = 'hcli';
    hcli{i}.parameter = bt.parameter;
    hcli{i}.parameter(options.free_pars) = p;
    hcli{i}.mesh = uniformmesh;
    hcli{i}.degree = options.degree;
    hcli{i}.profile = ups;
    hcli{i}.period = 2*T;
    hcli{i}.x1 = x0; % saddle
    hcli{i}.x2 = x0; % (homoclinic)

    %% lambda_{v,w}, v, w
    stst.kind='stst';
    stst.parameter=hcli{i}.parameter;
    stst.x=hcli{i}.x1;
    m=df_mthod('stst','cheb');
    stst.stability=p_stabil(funcs,stst,m.stability);
    if isempty(stst.stability.l1) || max(real(stst.stability.l1))<0
        error('P_TOHCLI: no unstable eigenmodes found');
    end
    lambda_v=stst.stability.l1(:);
    v_sel=real(lambda_v)>0;
    lambda_v=lambda_v(v_sel);
    hcli{i}.lambda_v=lambda_v;
    hcli{i}.v=stst.stability.v(:,v_sel);
    hcli{i}.lambda_w = hcli{i}.lambda_v;
    hcli{i}.w = hcli{i}.v;
    hcli{i}.w = hcli{i}.w/norm(hcli{i}.w); % normalize w to for the BVP

    %% alpha, espilon
    hcli{i}.alpha=hcli{i}.v\(hcli{i}.profile(:,1)-hcli{i}.x1);
    hcli{i}.epsilon=norm(hcli{i}.alpha);
    hcli{i}.alpha=hcli{i}.alpha/hcli{i}.epsilon;
end
