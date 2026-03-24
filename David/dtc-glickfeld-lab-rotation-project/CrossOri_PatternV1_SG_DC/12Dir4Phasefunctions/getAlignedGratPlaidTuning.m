% Input is avg_resp_dir, nCells x nDir x nMaskPhase x (1: grating, 2:
% plaid) x (1: mean resp, 2: std)

% This function aligns all grating and plaid responses to preferred grating
% direction (computed simply as max avg response across nDirs)

function [gratAligned, plaidAligned] = getAlignedGratPlaidTuning(avg_resp_dir, aligntype, plaid_offset_deg)

    if nargin < 2
        aligntype = '';
    end
    if nargin < 3 || isempty(plaid_offset_deg)
        plaid_offset_deg = 45;
    end
    %this allows us to specify that if we are aligning intracellular data,
    %we do not need to shift the plaid resp by 60 degrees, since these
    %responses are defined by their coherent direction.
    
    nCells      = size(avg_resp_dir,1);
    nMaskPhas   = size(avg_resp_dir,3);
   
    int = 360/size(avg_resp_dir,2);   % Direction step size, in degrees

    respgrat = avg_resp_dir(:,:,1,1,1);     % Get grating responses only
    [prefresp, prefdir] = max(respgrat,[],2);   % Find preferred grating direction
    vecshift = zeros(size(respgrat,1),1)+size(avg_resp_dir,2)/2;
    vecshift = vecshift - prefdir;
    respplaid = avg_resp_dir(:,:,:,2,1);
    
    gratAligned =[];
    plaidAligned = [];
    
    for iCell = 1:size(respgrat,1)
        gratAligned(iCell,:) = circshift(respgrat(iCell,:),vecshift(iCell),2);
        
        if strcmp(aligntype, 'whole-cell')
            for i = 1:nMaskPhas
                plaidAligned(iCell,:,i) = circshift(respplaid(iCell,:,i),vecshift(iCell),2); %align, but don't shift by 60 degrees
            end

        else
            for i = 1:nMaskPhas
                plaidAligned(iCell,:,i) = circshift(respplaid(iCell,:,i),vecshift(iCell)+((plaid_offset_deg/2)./int),2); %align, and shift by half plaid angle
            end
        end
    end

end