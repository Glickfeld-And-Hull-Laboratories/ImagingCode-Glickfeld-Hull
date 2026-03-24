% Input is avg_resp_dir, nCells x nDir x nMaskPhase x (1: grating, 2:
% plaid) x (1: mean resp, 2: std)

function [gratingFitStruct] = getGratingTuningCurveFit(avg_resp_dir)
    nCells      = size(avg_resp_dir,1);
    nStimDir    = size(avg_resp_dir,2);
    stimDirs    = 0:(360/nStimDir):359;
    stimDirs    = stimDirs;

    y_fits = zeros(360,nCells);
    dirs = deg2rad(0:1:359);
    for iCell = 1:nCells
        [b_hat_all(iCell,1), k1_hat_all(iCell,1), R1_hat_all(iCell,1), R2_hat_all(iCell,1), u1_hat_all(iCell,1), u2_hat_all(iCell,1), sse_all(iCell,1),R_square_all(iCell,1)] = miaovonmisesfit_dir(deg2rad(stimDirs),avg_resp_dir(iCell,:,1,1,1));
        dir_yfit_all(:,iCell) = b_hat_all(iCell,1)+R1_hat_all(iCell,1).*exp(k1_hat_all(iCell,1).*(cos(dirs-u1_hat_all(iCell,1))-1))+R2_hat_all(iCell,1).*exp(k1_hat_all(iCell,1).*(cos(dirs-u1_hat_all(iCell,1)-pi)-1));
    end

    gratingFitStruct.b      = b_hat_all;
    gratingFitStruct.k1     = k1_hat_all;
    gratingFitStruct.R1     = R1_hat_all;
    gratingFitStruct.R2     = R2_hat_all;
    gratingFitStruct.u1     = u1_hat_all;
    gratingFitStruct.u2     = u2_hat_all;
    gratingFitStruct.sse    = sse_all;
    gratingFitStruct.Rsq    = R_square_all;
    gratingFitStruct.yfit   = dir_yfit_all;
end