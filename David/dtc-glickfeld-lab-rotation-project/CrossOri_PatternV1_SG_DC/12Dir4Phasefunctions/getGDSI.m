function [gDSIstruct] = getGDSI(responses, directions_deg)
% getGDSI  Global direction selectivity index (Li et al. 2025, Curr Biol)
%
%   gDSIstruct = getGDSI(responses, directions_deg)
%
% Input:
%   responses      - nCells x nDir (mean response at each direction)
%   directions_deg - 1 x nDir (stimulus directions in degrees)
%
% Output:
%   gDSIstruct.gDSI        - 1 x nCells (range 0-1)
%   gDSIstruct.prefDir_deg - 1 x nCells (range 0-360, via atan2)
%   gDSIstruct.DS_ind      - indices where gDSI > 0.2

    nCells = size(responses, 1);
    dirs_rad = deg2rad(directions_deg(:)');

    gDSI = zeros(1, nCells);
    prefDir_deg = zeros(1, nCells);

    for iCell = 1:nCells
        R = responses(iCell, :);
        vec_sum_x = sum(R .* cos(dirs_rad));
        vec_sum_y = sum(R .* sin(dirs_rad));
        gDSI(iCell) = sqrt(vec_sum_x^2 + vec_sum_y^2) / sum(R);
        prefDir_deg(iCell) = mod(rad2deg(atan2(vec_sum_y, vec_sum_x)), 360);
    end

    gDSIstruct.gDSI = gDSI;
    gDSIstruct.prefDir_deg = prefDir_deg;
    gDSIstruct.DS_ind = find(gDSI > 0.2);
end
