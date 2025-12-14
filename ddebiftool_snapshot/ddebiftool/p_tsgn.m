function [delay_nr,tz]=p_tsgn(funcs,point)
%% find negative delays
% function [delay_nr,tz]=p_tsgn(point)
% INPUT:
%       funcs problem function
%       point a point 
% OUTPUT:
%	delay_nr number of negative delay 
%       tz (for psol only) we want tau(tz)=0 and dtau/dt(tz)=0
%
% (c) DDE-BIFTOOL v. 2.00, 30/11/2001
%
% $Id: p_tsgn.m 373 2019-09-03 22:13:51Z jansieber $
%
%%
d=dde_num_delays(funcs); % number of delays
delay_nr=[];
tz=[];
tau_all=p_tau(funcs,point,1:d);
switch point.kind
    case {'stst','hopf','fold'}
        delay_nr=find(tau_all<0,1,'first');
    case {'psol','coll'}
        %% compute roots of derivative of interpolation polynomial
        for i=1:d
            [t_ext,tauval]=dde_coll_roots(setfield(point,'profile',tau_all(i,:,1)),1,'diff',1); %#ok<SFLD>
            [min_delay,i_min]=min(tauval);
            if min_delay<=0
                tz=t_ext(i_min);
                delay_nr=i;
                break
            end
        end
    case 'hcli'
        error('P_TSGN: this routine is not (yet) implemented for connecting orbits.');
    otherwise
        error('P_TSGN: point kind %s not recognized.',point.kind);
end
end
