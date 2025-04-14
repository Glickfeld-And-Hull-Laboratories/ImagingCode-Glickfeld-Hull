%

% expected data = resp_cell; {nStimDir x nPhase x (1 for gratings or 2 for
% plaids)}(nCells x nTrials)

function [boot_base, boot_amp, boot_rsq, boot_sse, boot_Zp, boot_Zc] = bootstrap_onephase(data, nPhases, pd, nboots)

nStimDir = size(data,1);
nCells = size(data{1,1,1},1);

boot_base = [];
boot_amp = [];
boot_rsq = [];
boot_sse = [];
boot_Zc = zeros(nPhases,nCells,nboots);
boot_Zp = zeros(nPhases,nCells,nboots);

for i = 1:nboots
    
    resps = double.empty(nCells,0);
    avg_resp_new = zeros(nCells, nStimDir, nPhases);
    avg_resp_grat = zeros(nCells, nStimDir);

    for id = 1:nStimDir
        resps = data{id,1,2};
        for ip = 1:nPhases
            ntrials = round(random(pd));
            if ntrials > size(resps,2)
                ntrials = size(resps,2);
            end
            resp_cell_new{id,ip} = resps(:,randperm(size(resps,2),ntrials));
            avg_resp_new(:,id,ip) =  mean(resp_cell_new{id,ip}(:,:),2); 
        end
        avg_resp_grat(:,id) =  mean(data{id,1,1}(:,:),2);
    end

    component = avg_resp_grat(:,:)+circshift(avg_resp_grat(:,:),-120./nStimDir,2);
    pattern = circshift(avg_resp_grat(:,:),-60./nStimDir,2);

    %Calculate pattern and component prediction
    comp_corr = zeros(nPhases,nCells);
    patt_corr = zeros(nPhases,nCells);
    comp_patt_corr = zeros(nPhases,nCells);
    plaid_corr = zeros(1,nCells);
    plaid_corr_rand = zeros(1,nCells);

    for iCell = 1:nCells
        for ip = 1:nPhases
            comp_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,ip),component(iCell,:)));
            patt_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,ip),pattern(iCell,:)));
            comp_patt_corr(ip,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
        end
    end
    Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
    Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
    Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
    Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));

    PCI = Zp-Zc;

    phase = 0:(360/nPhases):359;
    phase_range = 0:359;
    for iCell = 1:nCells
        [b_hat_all(iCell,1), amp_hat_all(iCell,1), per_hat_all(iCell,1),pha_hat_all(iCell,1),sse_all(iCell,1),R_square_all(iCell,1)] = sinefit_PCI(deg2rad(phase),PCI(:,iCell));
        yfit_all(iCell,:,1) = b_hat_all(iCell,1)+amp_hat_all(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(iCell,1) + 2.*pi/pha_hat_all(iCell,1)));
    end
    boot_base = [boot_base, b_hat_all];
    boot_amp = [boot_amp, amp_hat_all];
    boot_rsq = [boot_rsq, R_square_all];
    boot_sse = [boot_sse, sse_all];
    boot_Zp(:,:,i)=Zp;
    boot_Zc(:,:,i)=Zc;
end

end


