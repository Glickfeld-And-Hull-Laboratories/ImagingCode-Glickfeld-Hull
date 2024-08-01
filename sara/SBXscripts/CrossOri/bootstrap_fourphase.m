%Four phase bootstrap -- shuffle trials across plaid phase

% expected data = resp_cell; {nStimDir x nPhase x (1 for gratings or 2 for
% plaids)}(nCells x nTrials)

function [boot_base, boot_amp, boot_Zc, boot_Zp] = bootstrap_fourphase(data, nboots)

nStimDir = size(data,1);
nMaskPhas = size(data,2);
nCells = size(data{1,1,1},1);

boot_base = [];
boot_amp = [];
boot_Zc = zeros(nMaskPhas,nCells,nboots);
boot_Zp = zeros(nMaskPhas,nCells,nboots);

for i = 1:nboots
	ntrial = {};
    resps = double.empty(nCells,0);
    avg_resp_dir_shuf = zeros(nCells, nStimDir, nMaskPhas);
    avg_resp_grat = zeros(nCells, nStimDir);

    for id = 1:nStimDir
        for ip = 1:nMaskPhas
            ntrial{id,ip} = size(data{id,ip,2},2);
            resps = [resps, data{id,ip,2}];
        end
        resps_shuf = resps(:,randperm(size(resps,2)));
        start = 1;
        t=0;
        for ip = 1:nMaskPhas
            t = t + ntrial{id,ip};
            resp_cell_shuf{id,ip} = resps_shuf(:,start:t);
            avg_resp_dir_shuf(:,id,ip) =  mean(resp_cell_shuf{id,ip}(:,:),2);
        end
        avg_resp_grat(:,id) =  mean(data{id,1,1}(:,:),2);
    end

    component = avg_resp_grat(:,:)+circshift(avg_resp_grat(:,:),-120./nStimDir,2);
    pattern = circshift(avg_resp_grat(:,:),-60./nStimDir,2);

    comp_corr = zeros(nMaskPhas,nCells);
    patt_corr = zeros(nMaskPhas,nCells);
    comp_patt_corr = zeros(nMaskPhas,nCells);

    for iCell = 1:nCells
        for ip = 1:nMaskPhas
            comp_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir_shuf(iCell,:,ip),component(iCell,:)));
            patt_corr(ip,iCell) = triu2vec(corrcoef(avg_resp_dir_shuf(iCell,:,ip),pattern(iCell,:)));
            comp_patt_corr(ip,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
        end
    end
    Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
    Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
    Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
    Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));

    PCI = Zp-Zc;

    phase = 0:(360/nMaskPhas):359;
    phase_range = 0:359;
    for iCell = 1:nCells
        [b_hat_all(iCell,1), amp_hat_all(iCell,1), per_hat_all(iCell,1),pha_hat_all(iCell,1),sse_all(iCell,1),R_square_all(iCell,1)] = sinefit_PCI(deg2rad(phase),PCI(:,iCell));
        yfit_all(iCell,:,1) = b_hat_all(iCell,1)+amp_hat_all(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(iCell,1) + 2.*pi/pha_hat_all(iCell,1)));
    end
    boot_base = [boot_base, b_hat_all];
    boot_amp = [boot_amp, amp_hat_all];
    boot_Zp(:,:,i)=Zp;
    boot_Zc(:,:,i)=Zc;
end
end


