%% test sco_
clear
format compact
syms x xp p 
brat=[xp; -p*exp(x)];
Fs=sco_sym2funcs(brat,{[x;xp],p},{'x','p'},'vector',[1,1],'maxorder',4);
%%
bratu={... % function
    @(x,p)[...
    x(2,:);...
    -p(1,:).*exp(x(1,:))],...
    ... % 1st derivative
    @(x,p,dx,dp)[...
    dx(2,:);...
    -dp(1,:).*exp(x(1,:))-p(1,:).*exp(x(1,:)).*dx(1,:)],...
    ... % 2nd derivative
    @(x,p,dx,dp)[...
    0*x(2,:);...
    -2*dp(1,:).*exp(x(1,:)).*dx(1,:)-p(1,:).*exp(x(1,:)).*dx(1,:).^2]};
F=sco_fun(bratu,{'x','p'});
F0=sco_fun(bratu{1},{'x','p'});
F0a=sco_fun(bratu(1),{'x','p'});
F1=sco_fun(bratu(1:2),{'x','p'});
funs={Fs,F,F0,F0a,F1};
n_funs=length(funs);
f_tests={'','x','p',{'x','x'},{'p','x'},{'x','p'},{'p','p'},0,1,2,...
    'x*v',{'x','x*v'}};
n_tests=length(f_tests);
nargs=2*ones(1,n_tests);
nargs(end-1:end)=3;
fun_comb=cell(n_funs,n_tests);
for i=1:n_tests
    for k=1:n_funs
        fun_comb{k,i}=funs{k}(f_tests{i});
    end
end
nvec=5;
rng(0);
x0=rand(2,nvec);
p0=rand(1);
v0=rand(2,nvec);
args={x0,p0,v0};
err=NaN(n_funs-1,n_tests);
for i=1:n_tests
    fprintf('==\n');
    disp(f_tests{i});
    for k=1:n_funs
        y{k}=fun_comb{k,i}(args{1:nargs(i)});
    end
    disp(cellfun(@(c){size(c)},y));
    for k=2:n_funs
        err(k-1,i)=norm(y{1}-y{k},'fro');
    end
    fprintf('err=%g\n',err(k-1,i));
end
