% Input is avg_resp_dir, nCells x nDir x nMaskPhase x (1: grating, 2:
% plaid) x (1: mean resp, 2: std)

% This function aligns all grating responses' SEM to preferred gratig
% direction (computed simply as max avg response across nDirs)

% Same as getAlignedGratPlaidTuning.m function

function [semGratAligned] = getAlignedGratSEM(avg_resp_dir)

    nCells      = size(avg_resp_dir,1);
    nMaskPhas   = size(avg_resp_dir,3);
   
    int = 360/size(avg_resp_dir,2);   % Direction step size, in degrees

    respgrat = avg_resp_dir(:,:,1,1,1);     % Get grating responses only
    gratsem = avg_resp_dir(:,:,1,1,2);
    [prefresp, prefdir] = max(respgrat,[],2);   % Find preferred grating direction
    vecshift = zeros(size(respgrat,1),1)+6; 
    vecshift = vecshift - prefdir;
    
    semGratAligned = [];
    
    for iCell = 1:size(respgrat,1)
        semGratAligned(iCell,:) = circshift(gratsem(iCell,:),vecshift(iCell),2);
    end

end