
% expected data = resp_cell; {nStimDir x nPhase x (1 for gratings or 2 for
% plaids)}(nCells x nTrials)

%ntrials -- expected to be array of list of ntrials (e.g., [12 24 36 48])

function [std_avg] = bootstrap_trials(data, num_trials, nboots)

nTrials = length(num_trials);

nStimDir = size(data,1);
nCells = size(data{1,1,1},1);

boot_PCI = zeros(nCells, nTrials, nboots);


for i = 1:nboots
    
    resps = double.empty(nCells,0);
    avg_resp_new = zeros(nCells, nStimDir, length(nTrials));
    avg_resp_grat = zeros(nCells, nStimDir, length(nTrials));
    
    for id = 1:nStimDir
        resps = data{id,1,2};
        for it = 1:nTrials    %:nPhases
            t = num_trials(it);
            resp_cell_new{id,it} = resps(:,randperm(size(resps,2),t));
            avg_resp_new(:,id,it) =  mean(resp_cell_new{id,it}(:,:),2);
        end
        avg_resp_grat(:,id) =  mean(data{id,1,1}(:,:),2);
    end
    
    component = avg_resp_grat(:,:)+circshift(avg_resp_grat(:,:),-120./nStimDir,2);
    pattern = circshift(avg_resp_grat(:,:),-60./nStimDir,2);

    %Calculate pattern and component prediction
    comp_corr = zeros(nTrials,nCells);
    patt_corr = zeros(nTrials,nCells);
    comp_patt_corr = zeros(nTrials,nCells);
    plaid_corr = zeros(1,nCells);
    plaid_corr_rand = zeros(1,nCells);
   
    for iCell = 1:nCells
        for it = 1:nTrials
            comp_corr(it,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,it),component(iCell,:)));
            patt_corr(it,iCell) = triu2vec(corrcoef(avg_resp_new(iCell,:,it),pattern(iCell,:)));
            comp_patt_corr(it,iCell) = triu2vec(corrcoef(component(iCell,:),pattern(iCell,:)));
        end
    end
    Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
    Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
    Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nStimDir-3));
    Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nStimDir-3));
    
    PCI = Zp-Zc;
    
    boot_PCI(:,:,i) = PCI';
end

std_allCells = std(boot_PCI,0,3);
std_avg = mean(std_allCells,1);

end




