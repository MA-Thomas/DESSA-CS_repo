function [wait_time,integPropValue] = computeWaitTime_preComputed_IntegratedPropensity(...
    IntegratedPropensity_curve, times_list, t_offset_b, Pk)

% This file assumes there is no nonzero t_offset_a for particle A.

% Note:  IntegratedPropensity_curve = rate_constant * CDF_F;


% TODO: 
%{
1. determine time at which IntegatedPropensity(time|d) == Pk + IntegatedPropensity(t_offset|d)
2. determine the wait_time (depends on whether the t_offset is non zero)

%}

if t_offset_b > 0
    [~,ind_t_offset] = min(abs(times_list - t_offset_b));
    IntegratedPropensity_at_t_offset = IntegratedPropensity_curve(ind_t_offset);
else
    IntegratedPropensity_at_t_offset = 0;
end

if IntegratedPropensity_curve(end) - IntegratedPropensity_at_t_offset < Pk
    % No reaction in allowed max diffusion time
    wait_time = -1;
    integPropValue = -1;
    return
end

inds = find(IntegratedPropensity_curve >= Pk + IntegratedPropensity_at_t_offset);



%if ~isempty(inds)
    
    wait_time = times_list(inds(1));
    integPropValue = IntegratedPropensity_curve( inds(1) );    

    assert(wait_time>0)
% else
%     % No reaction in allowed max diffusion time
%     wait_time = -1;
%     integPropValue = -1;
% 
% end

% % Do not consider diffusion times greater than <maxAllowed_t_wait>
% % which is set by the distance to the nearest relevant boundary.
% if wait_time > maxAllowed_t_wait; wait_time = -1; end









end

