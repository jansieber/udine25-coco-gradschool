function pat=dde_funformat_create(varargin)
default={'f',NaN,'x',NaN,'rhstau_args',NaN,'tauprovided',NaN,'par',NaN,...
    'funcs',[],'info',[]};%,...
    %'Jfpar',[],'Jfx',[],'Jtaux',[],'Jtaupar',[]};
pat=dde_set_options(default,varargin,'pass_on');
fn={'f','x','rhstau_args','tauprovided','par'};
for i=1:length(fn)
    pat.(fn{i})=reshape(pat.(fn{i}),1,[]);
end
end
