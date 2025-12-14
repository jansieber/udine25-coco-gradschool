function result_array=extract_from_tr(trPO_array,component,biftype,ip)
%% extract components from extended solution branch or point array
%
% $Id: extract_from_tr.m 309 2018-10-28 19:02:42Z jansieber $
%
%% check if input is branch rather than point array
if ~isfield(trPO_array,'kind') && isfield(trPO_array,'point')
    trPO_array=trPO_array.point;
end
dim=ip.dim;
type={'kind','solution','eigenvector','omega','solution_for_stability','vpoint'};
for i=1:length(trPO_array)
    trPO=trPO_array(i);
    switch component
        case type{1} % kind
            result_array=biftype;
            break
        case type{2} % solution
            result=trPO;
            result.profile=result.profile(1:dim,:);
            result.parameter=result.parameter(1:ip.nuserpar);
        case type{3} % eigenvector
            result=trPO;
            result.profile=result.profile(dim+1:end,:);
            result.parameter=zeros(1,ip.nuserpar);
            result.parameter(ip.nullparind(:,1))=result.parameter(ip.nullparind(:,2));
        case type{4} % omega
            result=trPO.parameter(ip.omega);
        case type{5} % including omega
            result=trPO;
            result.profile=result.profile(1:dim,:);
            result.parameter=result.parameter(1:ip.nuserpar);
            result.omega=trPO.parameter(end-1);
            result.flag=biftype;
        case type{6} % vpoint for initialization of nullvector
            result=trPO;
            result.profile=result.profile(dim+(1:ip.extdim),:);
            result.parameter=result.parameter(ip.nullparind(:,2).');
            result.period=trPO.parameter(ip.periodnullpar);
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