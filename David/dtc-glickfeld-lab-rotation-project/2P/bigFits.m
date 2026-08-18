%% big fits code

% Input is avg_resp_dir, nCells x nDir x nMaskPhase x (1: grating, 2:
% plaid) x (1: mean resp, 2: std)
function [DSIstruct, ZpZcStruct, plaid_corr, gratingFitStruct, ZpZcPWdist, phaseModStruct] = bigFits(avg_resp_dir)
    data_type = 'alignedTestDir';

    DSIstruct           = getDSIstruct(avg_resp_dir);
    ZpZcStruct          = getZpZcStruct(avg_resp_dir, data_type);
    plaid_corr          = getPlaidTuningCorrelations(avg_resp_dir);
    gratingFitStruct    = getGratingTuningCurveFit(avg_resp_dir);
    ZpZcPWdist          = getZpZcPWdist(ZpZcStruct);
    phaseModStruct      = get4PhaseModulationFit(ZpZcStruct);

end






