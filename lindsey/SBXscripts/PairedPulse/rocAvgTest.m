roc_avgCells = zeros(noff, nDelta,nexp);
for iexp = 1:nexp
     mouse = mouse_mat(iexp,:);
     date = date_mat(iexp,:);
     load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_roc180v23.mat']))
     load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
     load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_deltaResp.mat']))

    [x y] = sort(theta_90);
    tot = sum(~isnan(x));
    ind = y(1:tot);
    
    ppResp_avg = cell(noff+1,nDelta);
    for ioff= 1:noff+1
        for iDelta = 1:nDelta
            ppResp_avg{ioff,iDelta} = mean(ppResp{ioff,iDelta}(ind,:),1);
        end
    end


    for ioff = 1:noff
        for iDelta = 1:nDelta
            roc_avgCells(ioff,iDelta,iexp) = roc_gh(ppResp_avg{1,end}, ppResp_avg{ioff,iDelta});
        end
    end
end

figure; errorbar(mean(roc_avgCells,3)', std(roc_avgCells,[],3)'./sqrt(nexp))


plot(deltas, roc_avgCells(1,:)); hold on; plot(deltas, roc_avgCells(2,:));