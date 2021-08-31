function [c_50] = find_c50(con_resp_curve, cons)
%find_c50 
%   calculates the interpolated contrast at which the response is half max
%   enter a single vector of contrasts and responses
max_pref_resp = max(con_resp_curve);
half_max = max_pref_resp / 2;
c_50 = interp1(con_resp_curve, cons,half_max, 'linear');
end

