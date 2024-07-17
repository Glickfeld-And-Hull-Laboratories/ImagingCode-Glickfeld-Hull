%data_dfof is [nCells nTrials] 
%targetDelta is the stimulus condition on each trial (deltas is unique
%conditions, and nDelta is the number of conditions)
%nboot is typically 1000
for iC = 1:nCells
    fprintf('%d\n', iC)
    %first set up the average response to all trials, then bootstrap trials
    %with replacement
    for iboot = 1:nboot+1
        for idelta = 1:nDelta
            ind = intersect(find(targetDelta == deltas(idelta)));
            if iboot>1
                ind_use = ind(randsample(1:length(ind),length(ind),1));
            else
                ind_use = ind;
            end
            delta_resp(iC,idelta,iboot) = nanmean(data_dfof(iC,ind_use),2);
        end
    end
    %first fit tuning curve for average of all trials
    data = squeeze(delta_resp(iC,:,1));
    [b_hat(iC,1), k_hat(iC,1), R_hat(iC,1),u_hat(iC,1),sse(iC,1),R_square(iC,1)] = miaovonmisesfit_ori(deg2rad(theta),data);
    y_fit(:,iC,1) = b_hat+R_hat.*exp(k_hat(iC,1).*(cos(2.*(deg2rad(theta_smooth)-u_hat))-1));
    %measure the preferred orientation from the fit
    [y_max max_ori(iC,1,1)] = max(y_fit(:,iC,1),[],1);
    %fit tuning curve for all bootstraps
    for iboot = 2:nboot+1
        fprintf('.')
        data = squeeze(delta_resp(iC,:,iboot));
        [b_hat, k_hat_boot, R_hat,u_hat,sse,R_square(iC,iboot)] = miaovonmisesfit_ori(deg2rad(theta),data);
        y_fit(:,iC,iboot) = b_hat+R_hat.*exp(k_hat_boot.*(cos(2.*(deg2rad(theta_smooth)-u_hat))-1));
        [y_max max_ori(iC,1,iboot)] = max(y_fit(:,iC,iboot),[],1);
    end
    %plot average and fit for each cell
    subplot(n, n2, iCell)
    scatter(deg2rad(theta), squeeze(delta_resp(iC,:,1)),'ok');
    hold on; plot(deg2rad(0:1:180), y_fit(:,iC,1),'-r')
    title(num2str(chop(R_square(iC,1),2)))
end

theta_90 = nan(1,nCells); %this will be a metric of reliability- the larger it is, the less reliable
for iCell = 1:nCells
    if ~isnan(R_square(iC,1))
        %find the distance for all bootrstraps from the average
        theta_dist = abs(theta_smooth(squeeze(max_ori(iC,1)))-theta_smooth(squeeze(max_ori(iC,2:nboot+1))));
        %correct for circularity of ori (not needed for SF)
        theta_dist(find(theta_dist>90)) = 180-theta_dist(find(theta_dist>90));
        %sort distances to find the 90th percentile
        theta_sort = sort(theta_dist,'ascend');
        theta_90(1,iC) = theta_sort(ceil(nboot*.9));
    end
end

