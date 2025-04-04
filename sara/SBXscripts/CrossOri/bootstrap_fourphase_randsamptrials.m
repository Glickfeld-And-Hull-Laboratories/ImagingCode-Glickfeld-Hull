%Four phase bootstrap -- randomly sample trials for each stimulus
%condition with replacement

% expected data = resp_cell; {nStimDir x nPhase x (1 for gratings or 2 for
% plaids)}(nCells x nTrials)


function [boot_maxInd_all] = bootstrap_fourphase_randsamptrials(data, nboots)

nStimDir = size(data,1);
nMaskPhas = size(data,2);
nCells = size(data{1,1,1},1);

boot_maxInd_all = [];

for i = 1:nboots
	ntrial = {};
    resps = double.empty(nCells,0);
    avg_resp_dir_shuf = zeros(nCells, nStimDir, nMaskPhas);
    avg_resp_grat = zeros(nCells, nStimDir);

    for id = 1:nStimDir
        ntrial = size(data{id,1,1},2);
        resps = data{id,1,1};
        r_trials = randsample(ntrial,ntrial,true);
        r_resps = resps(:,r_trials);
        r_resp_cell{id} = r_resps;
        r_avg_resp_dir(id,:) = mean(r_resps,2);
    end

    [r_grat_dir_max r_grat_dir_maxInd] = max(r_avg_resp_dir,[],1);
    boot_maxInd_all = [boot_maxInd_all, r_grat_dir_maxInd'];
end
end


