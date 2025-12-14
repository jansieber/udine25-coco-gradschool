%% check if parameter bound has been crossed
function varargout=dde_check_parameterbounds(funcs,mth,pts,free_par,min_bound,max_bound,output)
dostop=false;
bdlist={min_bound,max_bound};
sg=[-1,1];
bd=struct('ind',NaN,'frac',NaN,'param',NaN);
for i=1:length(bdlist)
    bds=bdlist{i};
    for j=1:size(bds,1)
        if sg(i)*(pts(end).parameter(bds(j,1))-bds(j,2))>0 % over maximum
            dostop=true;
            iparhit=j;
            iminmaxhit=i;
            if isfield(mth.continuation,'warnings') && ...
                    mth.continuation.warnings && strcmp(output,'flag')
                fprintf('parameter bound %g for parameter(%d) reached\n',...
                    bds(iparhit,2),bds(iparhit,1));
            end
            break
        end
    end
    if dostop
        break
    end
end
if strcmp(output,'flag')
    varargout={dostop};
    return
end
if ~dostop || length(pts)<2
    varargout={pts(end),false};
    return
end
%% construct initial guess
bd.ind=bds(iparhit,1);
pout=pts(end).parameter(bds(iparhit,1));
pin=pts(end-1).parameter(bds(iparhit,1));
bd.frac=(bds(iparhit,2)-pin)/(pout-pin);
bd.param=bds(iparhit,2);
secant=p_axpy(-1,pts(end-1),pts(end));
pts(end)=p_axpy(bd.frac,secant,pts(1));
%% correct final point
pini=pts(end);
bdlist{iminmaxhit}=bdlist{iminmaxhit}(iparhit,:);
bdlist{3-iminmaxhit}=repmat([0,0],0,1);
fcn=dde_funcs_add_cond(funcs,dde_cond_parameter('parameterbounds',bdlist{:}));
mth.point.extra_condition=true;
pref=pts(end-1);
[pout,suc]=p_correc(fcn,pini,free_par,[],mth.point,0,pref);
varargout={pout,suc};
end
