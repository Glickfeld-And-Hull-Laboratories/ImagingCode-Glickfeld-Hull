% Inputs expected to be 4 phases x nCells
% Output is 6 x nCells
function [ZpZcPWdist] = getZpZcPWdist(ZpZcStruct)
    Zp = ZpZcStruct.Zp;
    Zc = ZpZcStruct.Zc;
    nCells = size(Zp,2);
    ZpZc = [];

    for iCell = 1:nCells
        ZpZc(:,1,iCell) = Zp(:,iCell);
        ZpZc(:,2,iCell) = Zc(:,iCell);
    end
    
    ZpZcPWdist = double.empty(6,0);
    for iCell = 1:nCells
        ZpZcPWdist(:,iCell) = pdist(ZpZc(:,:,iCell));
    end
end