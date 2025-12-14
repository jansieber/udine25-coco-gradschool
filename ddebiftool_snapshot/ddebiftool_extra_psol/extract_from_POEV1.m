function result_array=extract_from_POEV1(poev1_array,component,ip)
%% extract components from POEV1 solution branch or point array
%
%
%% check if input is branch rather than point array
if ~isfield(poev1_array,'kind') && isfield(poev1_array,'point')
    poev1_array=poev1_array.point;
end
%% extract named components
dim=ip.dim;
npar=ip.nuserpar;
type={'kind','solution','nullvector','solution_for_stability',...
    'xtau_ind','eigenvector','vpoint'};
for i=1:length(poev1_array)
    pfold=poev1_array(i);
    switch component
        case type{1} %'kind'
            result_array='POfold';
            break
        case type{2} %'solution'
            result=pfold;
            result.profile=result.profile(1:dim,:);
            result.parameter=result.parameter(1:npar);
        case type([3,6]) %'nullvector','eigenvector'
            result=pfold;
            result.profile=result.profile(dim+(1:ip.extdim),:);
            result.period=result.parameter(ip.periodnullpar);
            result.parameter=zeros(1,npar);
            result.parameter(ip.nullparind(:,1))=result.parameter(ip.nullparind(:,2));
        case type{4} % including flag
            result=pfold;
            result.profile=result.profile(1:dim,:);
            result.parameter=result.parameter(1:npar);
            result.flag='POfold';
        case type{7} % vpoint for initialization of nullvector
            result=pfold;
            result.profile=result.profile(dim+(1:ip.extdim),:);
            result.parameter=result.parameter(ip.nullparind(:,2).');
            result.period=pfold.parameter(ip.periodnullpar);
        case type{5} % xtau_ind
            result_array=ip.ext_tau;
            break
        otherwise
            fprintf('known component types:\n');
            for k=1:length(type)
                fprintf('%s\n',type{k});
            end
            result_array=[];
            break;
    end
    result_array(i)=result; %#ok<AGROW>
end
end
