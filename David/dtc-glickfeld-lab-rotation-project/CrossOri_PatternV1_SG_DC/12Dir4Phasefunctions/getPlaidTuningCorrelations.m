% Input is avg_resp_dir, nCells x nDir x nMaskPhase x (1: grating, 2:
% plaid) x (1: mean resp, 2: std)

function    [plaid_corr]   = getPlaidTuningCorrelations(avg_resp_dir)
    nCells = size(avg_resp_dir,1);

    plaid_corr = zeros(1,nCells);
    plaid_corr1 = zeros(1,nCells);
    plaid_corr2 = zeros(1,nCells);
    plaid_corr3 = zeros(1,nCells);
    plaid_corr4 = zeros(1,nCells);
    plaid_corr5 = zeros(1,nCells);
    plaid_corr6 = zeros(1,nCells);
    
    for iCell = 1:nCells
        plaid_corr1(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,2,2,1)));
        plaid_corr2(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,3,2,1)));
        plaid_corr3(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,1,2,1),avg_resp_dir(iCell,:,4,2,1)));
        plaid_corr4(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,2,1),avg_resp_dir(iCell,:,3,2,1)));
        plaid_corr5(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,2,1),avg_resp_dir(iCell,:,4,2,1)));
        plaid_corr6(1,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,3,2,1),avg_resp_dir(iCell,:,4,2,1)));
        plaid_corr(1,iCell) = (plaid_corr1(1,iCell)+plaid_corr2(1,iCell)+plaid_corr3(1,iCell)+plaid_corr4(1,iCell)+plaid_corr5(1,iCell)+plaid_corr6(1,iCell))/6;
    end

end

