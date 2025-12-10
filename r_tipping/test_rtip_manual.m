%% test r_tipping_manual derivatives
pnames={'r','cm','beta','omega','cp','a','phi'}; % parameter names
[ip,npars]=structind_from_names(pnames);
rhs=@(ord,varargin)r_tipping_rhs_manual(ord,ip,varargin{:});
rhs_c={@(y,p)rhs(0,y,p),@(y,p,dy,dp)rhs(1,y,p,dy,dp),...
    @(y,p,dy,dp)rhs(2,y,p,dy,dp)};
Fs=sco_gen(@r_tipping_rhs);
F=sco_fun(rhs_c,{'x','par'});
F0=sco_fun(rhs_c{1},{'x','par'});
F0a=sco_fun(rhs_c(1),{'x','par'});
F1=sco_fun(rhs_c(1:2),{'x','par'});
funs={Fs,F,F0,F0a,F1};
n_funs=length(funs);
f_tests={'','x','par',{'x','x'},{'par','x'},{'x','par'},{'par','par'},0,1,2,...
    {'x*v'},{'x','x*v'},{'par','x*v'}};
n_tests=length(f_tests);
nargs=2*ones(1,n_tests);
nargs(end-2:end)=3;
fun_comb=cell(n_funs,n_tests);
for i=1:n_tests
    for k=1:n_funs
        fun_comb{k,i}=funs{k}(f_tests{i});
    end
end
nvec=1;
rng(0);
x0=rand(3,nvec);
p0=rand(length(pnames),1);
v0=rand(3,nvec);
args={x0,p0,v0};
err=NaN(n_funs-1,n_tests);
clear y
for i=1:n_tests
    fprintf('==\n');
    disp(f_tests{i});
    for k=1:n_funs
        y{k,i}=fun_comb{k,i}(args{1:nargs(i)});
    end
    disp(cellfun(@(c){size(c)},y(:,i)'));
    for k=2:n_funs
        err(k-1,i)=norm(y{1,i}-y{k,i},'fro');
    end
    fprintf('err=%g\n',err(k-1,i));
end
