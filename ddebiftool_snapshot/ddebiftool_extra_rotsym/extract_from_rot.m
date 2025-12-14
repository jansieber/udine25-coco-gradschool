function result_array=extract_from_rot(rot_array,component,dim)
%% extract components from rotating or modulated wave branch or point array
%
%
%% check if input is branch rather than point array
if ~isfield(rot_array,'kind') && isfield(rot_array,'point')
    rot_array=rot_array.point;
end
%% extract named components
type={'kind','solution','solution_for_stability','omega'};
for i=1:length(rot_array)
    rot=rot_array(i);
    if isfield(rot,'x')
       sol='x';
    else
        sol='profile';
    end
    switch component
        case type{1} %'kind'
            result_array='rot';
            break
        case type([2,3]) %'solution'
            result=rot;
            x=result.(sol);
            x=x(1:dim,:);
            result.(sol)=x;
            result.parameter=result.parameter(1:end-1);
        case type{4} % rotation speed
            result=rot.parameter(end);
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
