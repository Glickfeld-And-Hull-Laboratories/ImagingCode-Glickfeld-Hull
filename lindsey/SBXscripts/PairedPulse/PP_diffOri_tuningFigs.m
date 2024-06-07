close all
clear all
date = '240529';
ImgFolder = [{'002'}];
mouse = 'i1395';
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey\Analysis\2P';

load(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
load(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
load(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))

%% Find responsive cells
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey\Analysis\2P';

ind_base = setdiff(1:nTrials,ind_con);
[x1,y1] = ttest(squeeze(mean(data_dfof(base_win,:,1,ind_base),1))',squeeze(mean(data_dfof(resp_win,:,1,ind_base),1))','tail','left','alpha',0.05);
[x2,y2] = ttest(squeeze(mean(data_dfof(base_win,:,2,:),1))',squeeze(mean(data_dfof(resp_win,:,2,:),1))','tail','left','alpha',0.05);

h1 = zeros(ndir,nCells);
p1 = zeros(ndir,nCells);
h2 = zeros(nDelta,nCells);
p2 = zeros(nDelta,nCells);

data_dfof_dir = nan(size(data_dfof,1), nCells,ndir);
data_dfof_delta = nan(size(data_dfof,1), nCells,nDelta);
if ndir>1
    for idir = 1:ndir
        ind = setdiff(find(baseDir == dirs(idir)),ind_con);
        data_dfof_dir(:,:,idir) = nanmean(data_dfof(:,:,1,ind),4);
        for iCell = 1:nCells
            [h1(idir,iCell), p1(idir,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,1,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,1,ind),1)),'tail','left','alpha',0.05./(ndir-1));
        end
    end
else
    data_dfof_dir(:,:,1) = nanmean(data_dfof(:,:,1,:),4);
end
if nDelta>1
    for idelta = 1:nDelta
        ind = find(targetDelta == deltas(idelta));
        data_dfof_delta(:,:,idelta) = nanmean(data_dfof(:,:,2,ind),4);
        for iCell = 1:nCells
            [h2(idelta,iCell), p2(idelta,iCell)] = ttest(squeeze(nanmean(data_dfof(base_win,iCell,2,ind),1)),squeeze(nanmean(data_dfof(resp_win,iCell,2,ind),1)),'tail','left','alpha',0.05./(nDelta-1));
        end
    end
else
    data_dfof_delta(:,:,1) = nanmean(data_dfof(:,:,2,:),4);
end
save(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respSig.mat']), 'p2','h2')

good_ind_temp = find(x1+x2+sum(h1,1)+sum(h2,1));

%find late responding cells and remove
tc_all = squeeze(nanmean(mean(bsxfun(@minus,data_dfof(:,:,:,:),mean(data_dfof(base_win,:,:,:),1)),3),4));
resp_diff = diff(tc_all);
[max_val max_time] = max(resp_diff(20:end,:),[],1);

ind1 = find(max_time<base_win(end)-20);
ind2 = find(max_time>resp_win(end)-20);
ind3 = setdiff(1:nCells, [ind1 ind2]);
good_ind = intersect(good_ind_temp,ind3);

figure;
subplot(2,2,1)
plot(tc_all(:,ind1))
subplot(2,2,2)
plot(tc_all(:,ind2))
subplot(2,2,3)
plot(tc_all(:,ind3))
subplot(2,2,4)
plot(tc_all(:,good_ind))


%% account for blank trials
if length(ind_con)>0
    off_all = [offs; unique(cell2mat(input.tItiWaitFrames))];
    noff_all = length(off_all);
else
    off_all = offs;
    noff_all = noff;
end

%count trials
ind_n = zeros(noff_all, nDelta,ndir);
for ioff = 1:noff_all
    if length(ind_con)== 0 || ioff<noff_all 
        ind = setdiff(find(tFramesOff == off_all(ioff)),ind_con);
    else
        ind = ind_con;
    end
    for idelta = 1:nDelta
        ind_delta = find(targetDelta == deltas(idelta));
        for idir = 1:ndir
            ind_n(ioff,idelta,idir) = length(intersect(intersect(ind,ind_delta),find(baseDir == dirs(idir))));
        end
    end
end

%% find fit tau and subtract base response
tt = 1-20:1+99;

data_dfof_sub = nan(size(data_dfof,1),nCells,nTrials);
data_dfof_sub(:,:,ind_con) = data_dfof(:,:,2,ind_con);
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    fprintf('%d\n',iC)
    if find(find(x1+sum(h1,1)) == iC)
        tc = nanmean(data_dfof(resp_win(1):end,iC,2,ind_con),4);
        [A, tau, R_square] = singExpDecayFit(tc);
        if isnan(A)
            continue
        end
        fit_resp = A.*exp((1:length(tc))./tau);
        fit_resp_fill = [zeros(1,size(data_dfof,1)-length(fit_resp)) fit_resp];
        for ioff = 1:noff
            delay = -offs(ioff)-3;
            fit_shift = circshift(fit_resp_fill',delay);
            ind = setdiff(find(tFramesOff == off_all(ioff)),ind_con);
            for idelta = 1:nDelta
                ind2 = intersect(ind,find(targetDelta == deltas(idelta)));
                for itrial = 1:length(ind2)
                    it = ind2(itrial);
                    base_win_amp = mean(data_dfof(resp_win,iC,1,it),1)-mean(data_dfof(base_win,iC,1,it),1);
                    fit_resp_scale = (base_win_amp./fit_resp(1)).*fit_shift;
                    data_dfof_sub(:,iC,it) = data_dfof(:,iC,2,it)-fit_resp_scale;
                end
            end
        end
    else
        data_dfof_sub(:,iC,setdiff(1:nTrials,ind_con)) = data_dfof(:,iC,2,setdiff(1:nTrials,ind_con));
    end
end

figure;
[n n2] = subplotn(length(good_ind));
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    subplot(n,n2,iCell)
    ind = setdiff(find(tFramesOff == off_all(1)),ind_con);
    plot(squeeze(nanmean(data_dfof(:,iC,2,ind),4)))
    hold on
    plot(squeeze(nanmean(data_dfof_sub(:,iC,ind),3)))
end

save(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_decaySub.mat']), 'data_dfof_sub', 'good_ind', 'off_all', 'noff_all')
%% measure tuning function and bootstrap
    
nboot = 1000;
theta_smooth = (0:1:180);
theta = deltas;
[n n2] = subplotn(length(good_ind));

delta_resp = nan(nCells,nDelta,noff_all,1+nboot);
delta_resp_sub = nan(nCells,nDelta,noff_all,1+nboot);
max_dir = nan(nCells,1+nboot);
max_dir_sub = nan(nCells,1+nboot);
k_hat = nan(nCells,noff_all);
k_hat_sub = nan(nCells,noff_all);
y_fit = nan(length(theta_smooth),nCells,noff_all,nboot+1);
y_fit_sub = nan(length(theta_smooth),nCells,noff_all,nboot+1);
max_ori = nan(nCells,noff_all,nboot+1);
max_ori_sub = nan(nCells,noff_all,nboot+1);
R_square = nan(nCells,noff_all,nboot+1);

figure;
for iCell = 1:length(good_ind)    
    iC = good_ind(iCell);
    fprintf('%d\n', iC)
    if ~isnan(data_dfof_sub(1,iC,1))
        for ioff = 1:noff_all
            if length(ind_con)== 0 || ioff<noff_all 
                ind = setdiff(find(tFramesOff == off_all(ioff)),ind_con);
            else
                ind = ind_con;
            end
%             if ioff == noff_all
                for iboot = 1:nboot+1
                    for idelta = 1:nDelta
                        ind2 = intersect(ind,find(targetDelta == deltas(idelta)));
                        if iboot>1
                            ind_use = ind2(randsample(1:length(ind2),length(ind2),1));
                        else
                            ind_use = ind2;
                        end
                        delta_resp(iC,idelta,ioff,iboot) = nanmean(mean(data_dfof(resp_win,iC,2,ind_use),1)-mean(data_dfof(base_win,iC,2,ind_use),1),4);
                        delta_resp_sub(iC,idelta,ioff,iboot) = nanmean(mean(data_dfof_sub(resp_win,iC,ind_use),1)-mean(data_dfof_sub(base_win,iC,ind_use),1),3);
                    end
                end
%             else
%                 for idelta = 1:nDelta
%                     ind2 = intersect(ind,find(targetDelta == deltas(idelta)));
%                     delta_resp(iC,idelta,ioff,1) = nanmean(mean(data_dfof(resp_win,iC,2,ind2),1)-mean(data_dfof(base_win,iC,2,ind2),1),4);
%                     delta_resp_sub(iC,idelta,ioff,1) = nanmean(mean(data_dfof_sub(resp_win,iC,ind2),1)-mean(data_dfof_sub(base_win,iC,ind2),1),3);
%                 end
%             end
            data = squeeze(delta_resp(iC,:,ioff,1));
            data_sub = squeeze(delta_resp_sub(iC,:,ioff,1));
            [b_hat, k_hat(iC,ioff), R_hat,u_hat,sse(iC,ioff),R_square(iC,ioff,1)] = miaovonmisesfit_ori(deg2rad(theta),data);
            y_fit(:,iC,ioff,1) = b_hat+R_hat.*exp(k_hat(iC,ioff).*(cos(2.*(deg2rad(theta_smooth)-u_hat))-1));
            [y_max max_ori(iC,ioff,1)] = max(y_fit(:,iC,ioff,1),[],1);
            [b_hat, k_hat_sub(iC,ioff), R_hat,u_hat,sse_sub(iC,ioff),R_square_sub(iC,ioff,1)] = miaovonmisesfit_ori(deg2rad(theta),data_sub);
            y_fit_sub(:,iC,ioff,1) = b_hat+R_hat.*exp(k_hat_sub(iC,ioff).*(cos(2.*(deg2rad(theta_smooth)-u_hat))-1));
            [y_max max_ori_sub(iC,ioff,1)] = max(y_fit_sub(:,iC,ioff,1),[],1);
%             if ioff == noff_all
                [max_resp max_dir(iC,1)] = max(delta_resp(iC,:,ioff,1),[],2);
                for iboot = 2:nboot+1
                    fprintf('.')
                    [max_resp max_dir(iC,iboot)] = max(delta_resp(iC,:,ioff,iboot),[],2);
                    data = squeeze(delta_resp(iC,:,ioff,iboot));
                    [b_hat, k_hat_boot, R_hat,u_hat,sse,R_square(iC,ioff,iboot)] = miaovonmisesfit_ori(deg2rad(theta),data);
                    y_fit(:,iC,ioff,iboot) = b_hat+R_hat.*exp(k_hat_boot.*(cos(2.*(deg2rad(theta_smooth)-u_hat))-1));
                    [y_max max_ori(iC,ioff,iboot)] = max(y_fit(:,iC,ioff,iboot),[],1);
                    data_sub = squeeze(delta_resp_sub(iC,:,ioff,iboot));
                    [b_hat, k_hat_boot, R_hat,u_hat,sse,R_square_sub(iC,ioff,iboot)] = miaovonmisesfit_ori(deg2rad(theta),data_sub);
                    y_fit_sub(:,iC,ioff,iboot) = b_hat+R_hat.*exp(k_hat_boot.*(cos(2.*(deg2rad(theta_smooth)-u_hat))-1));
                    [y_max max_ori_sub(iC,ioff,iboot)] = max(y_fit_sub(:,iC,ioff,iboot),[],1);
                end
                subplot(n, n2, iCell)
                scatter(deg2rad(theta), squeeze(delta_resp(iC,:,ioff,1)),'ok');
                hold on; plot(deg2rad(0:1:180), y_fit(:,iC,ioff,1),'-r')
                title(num2str(chop(R_square(iC,ioff,1),2)))
%             end
        end
    end
end 
suptitle([date ' ' mouse])
print(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_tuningFits.pdf']),'-dpdf','-bestfit')

theta_90 = nan(noff_all,nCells);
max_dir_n = nan(1,nCells);
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    if ~isnan(R_square(iC,noff_all,1))
        max_dir_n(1,iC) = length(find(max_dir(iC,2:end)==max_dir(iC,1)));
        for ioff = 1:noff_all
            theta_dist = abs(theta_smooth(squeeze(max_ori(iC,ioff,1)))-theta_smooth(squeeze(max_ori(iC,ioff,2:nboot+1))));
            theta_dist(find(theta_dist>90)) = 180-theta_dist(find(theta_dist>90));
            theta_sort = sort(theta_dist,'ascend');
            theta_90(ioff,iC) = theta_sort(ceil(nboot*.9));
        end
    end
end

OSI = nan(nCells,noff_all);
pref_ori = nan(nCells,noff_all);
OSI_sub = nan(nCells,noff_all);
pref_ori_sub = nan(nCells,noff_all);
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    if theta_90(3,iC) <= 22.5
        for ioff = 1:noff_all
            pref_ori(iC,ioff) = theta_smooth(max_ori(iC,ioff,1));
            null_ori = max_ori(iC,ioff,1)-90;
            if null_ori <= 0
                null_ori= 180+null_ori;
            end
            pref_val = y_fit(max_ori(iC,ioff,1),iC,ioff,1);
            null_val = y_fit(null_ori,iC,ioff,1);
            if null_val < 0
                null_val = 0;
            end
            OSI(iC,ioff) = (pref_val-null_val)./ (pref_val+null_val);
        end
        for ioff = 1:noff_all
            pref_ori_sub(iC,ioff) = theta_smooth(max_ori_sub(iC,ioff,1));
            null_ori = max_ori_sub(iC,ioff,1)-90;
            if null_ori <= 0
                null_ori= 180+null_ori;
            end
            pref_val = y_fit_sub(max_ori_sub(iC,ioff,1),iC,ioff,1);
            null_val = y_fit_sub(null_ori,iC,ioff,1);
            if null_val < 0
                null_val = 0;
            end
            OSI_sub(iC,ioff) = (pref_val-null_val)./ (pref_val+null_val);
        end
    end
end

OSI_k = 1-exp(-2.*k_hat);
OSI_k_sub = 1-exp(-2.*k_hat_sub);
HWHM = 0.5.*acos((log(0.5)+k_hat)./k_hat);
HWHM_sub = 0.5.*acos((log(0.5)+k_hat_sub)./k_hat_sub);

save(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_deltaResp.mat']), 'good_ind', 'delta_resp', 'theta_90', 'max_dir', 'max_dir_n', 'pref_ori', 'OSI', 'y_fit', 'delta_resp_sub', 'pref_ori_sub', 'OSI_sub', 'y_fit_sub', 'k_hat', 'k_hat_sub', 'OSI_k', 'OSI_k_sub', 'HWHM', 'HWHM_sub', 'R_square', 'R_square_sub','sse','sse_sub')

%% test roc analysis
trialInd = cell(noff_all, nDelta);
ppResp = cell(noff_all, nDelta);
ppResp_sub = cell(noff_all, nDelta);
for ioff = 1:noff_all
    if length(ind_con)== 0 || ioff<noff_all 
        ind = setdiff(find(tFramesOff == off_all(ioff)),ind_con);
    else
        ind = ind_con;
    end
    for idelta = 1:nDelta
        ind2 = intersect(ind,find(targetDelta == deltas(idelta)));
        trialInd{ioff,idelta} = ind2;
        ppResp{ioff, idelta} = squeeze(mean(data_dfof(resp_win,:,2,ind2),1)-mean(data_dfof(base_win,:,2,ind2),1));
        ppResp_sub{ioff, idelta} = squeeze(mean(data_dfof_sub(resp_win,:,ind2),1)-mean(data_dfof_sub(base_win,:,ind2),1));
    end
end

roc_resp = nan(nCells,noff,nDelta);
roc_resp_sub = nan(nCells,noff,nDelta);
noise = ppResp{1,8};
noise_sub = ppResp_sub{1,8};
for iCell = 1:length(good_ind)
    iC = good_ind(iCell);
    for ioff = 1:noff
        for idelta = 1:nDelta
            temp_sig = ppResp{ioff,idelta};
            roc_resp(iC,ioff,idelta) = roc_gh(noise(iC,:),temp_sig(iC,:));
            temp_sig = ppResp_sub{ioff,idelta};
            roc_resp_sub(iC,ioff,idelta) = roc_gh(noise_sub(iC,:),temp_sig(iC,:));
        end
    end
end

save(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_roc180v23.mat']), 'roc_resp', 'roc_resp_sub','ppResp','ppResp_sub','trialInd')

%% avg timecourses
data_dfof_maxdir = nan(120,nCells);
for iC = 1:nCells
    if ~isnan(max_dir(iC,1))
        ind_use = intersect(ind_con, find(targetDelta == deltas(max_dir(iC,1))));
        data_dfof_maxdir(:,iC) = nanmean(data_dfof(:,iC,2,ind_use),4);
    end
end
tt= (1-prewin_frames:postwin_frames).*(1000/double(frameRateHz));
figure;
subplot(2,2,1)
plot(tt,data_dfof_maxdir)
xlabel('Time (s)')
ylabel('dF/F')
title('Max dir')
subplot(2,2,2)
plot(tt,nanmean(data_dfof_maxdir,2))
xlabel('Time (s)')
ylabel('dF/F')

ind1 = setdiff(find(tFramesOff == off_all(1)),ind_con);
ind2 = find(targetDelta == deltas(end));
ind_use = intersect(ind1,ind2);
cell_use = intersect(good_ind,find(nanmean(nanmean(data_dfof(resp_win,:,1,ind_use),1),4)>0.03));
subplot(2,2,3)
plot(tt,nanmean(data_dfof(:,cell_use,1,ind_use),4)')
xlabel('Time (s)')
ylabel('dF/F')
title('Repeat stim')
subplot(2,2,4)
plot(tt,nanmean(nanmean(data_dfof(:,cell_use,1,ind_use),4),2)')
xlabel('Time (s)')
ylabel('dF/F')
sgtitle([mouse ' ' date ])
print(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_avgTCs.pdf']),'-dpdf','-bestfit')

%% cell by cell figures
for iCell = 1:5 %length(good_ind)
    iC = good_ind(iCell);
%     if find(OSI_ind == iC)
%         sig_str = 'sig';
%     else
%         sig_str = 'non-sig';
%     end
    figure;
    for ioff = 1:noff_all
        if length(ind_con)== 0 || ioff<noff_all 
            ind = setdiff(find(tFramesOff == off_all(ioff)),ind_con);
        else
            ind = ind_con;
        end
        for idelta = 1:nDelta
            ind2 = intersect(ind,find(targetDelta == deltas(idelta)));
            subplot(noff_all,2,((ioff-1)*2)+1)
            plot(tt, nanmean(bsxfun(@minus,data_dfof(:,iC,2,ind2),mean(data_dfof(base_win,iC,2,ind2),1)),4))
            hold on
            subplot(noff_all,2,((ioff-1)*2)+2)
            errorbar(deltas(idelta),nanmean(mean(data_dfof(resp_win,iC,2,ind2),1)-mean(data_dfof(base_win,iC,2,ind2),1),4),nanstd(mean(data_dfof(resp_win,iC,2,ind2),1)-mean(data_dfof(base_win,iC,2,ind2),1),[],4)./sqrt(length(ind2)),'ok')
            hold on
        end
        plot(theta_smooth, squeeze(y_fit_all(:,iC,ioff,1)))
        subplot(noff_all,2,((ioff-1)*2)+1)
        title(num2str(off_all(ioff)*frameRateHz))
%         if ioff == noff_all
%             subplot(noff_all,2,((ioff-1)*2)+2)
%             title([num2str(deltas(max_dir(iC,:))) '- ' sig_str])
%         end
    end
    suptitle([mouse ' ' date '- Cell #' num2str(iC)]) 
    print(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_testResp_byDelta_Cell' num2str(iC) '.pdf']),'-dpdf')
end


for iCell = 1:5%length(good_ind)
    iC = good_ind(iCell);
%     if find(OSI_ind == iC)
%         sig_str = 'sig';
%     else
%         sig_str = 'non-sig';
%     end
    figure;
    for ioff = 1:noff_all
        if length(ind_con)== 0 || ioff<noff_all 
            ind = setdiff(find(tFramesOff == off_all(ioff)),ind_con);
        else
            ind = ind_con;
        end
        data = delta_resp_trial{iC,ioff,1};
        scatter(data(:,2), data(:,1),'o')
        hold on
        for idelta = 1:nDelta
            ind2 = intersect(ind,find(targetDelta == deltas(idelta)));
            subplot(noff_all,2,((ioff-1)*2)+1)
            %plot(deltas(idelta), nanmedian(mean(data_dfof_sub(resp_win,iC,ind2),1)-mean(data_dfof_sub(base_win,iC,ind2),1),3),'ok')
            %errorbar(deltas(idelta),nanmean(mean(data_dfof_sub(resp_win,iC,ind2),1)-mean(data_dfof_sub(base_win,iC,ind2),1),3),nanstd(mean(data_dfof_sub(resp_win,iC,ind2),1)-mean(data_dfof_sub(base_win,iC,ind2),1),[],3)./sqrt(length(ind2)),'ok')
            hold on
            subplot(noff_all,2,((ioff-1)*2)+2)
            errorbar(deltas(idelta),nanmean(mean(data_dfof_sub(resp_win,iC,ind2),1)-mean(data_dfof_sub(base_win,iC,ind2),1),3),nanstd(mean(data_dfof_sub(resp_win,iC,ind2),1)-mean(data_dfof_sub(base_win,iC,ind2),1),[],3)./sqrt(length(ind2)),'ok')
            hold on
        end
        plot(theta_smooth, squeeze(y_fit_all(:,iC,ioff,1)))
        ylim([-0.05 0.4])
        subplot(noff_all,2,((ioff-1)*2)+1)
        plot(theta_smooth, squeeze(y_fit_trial(:,iC,ioff,1)))
        title(num2str(off_all(ioff)*frameRateHz))
        ylim([-0.05 0.4])
%         if ioff == noff_all
%             subplot(noff_all,2,((ioff-1)*2)+2)
%             title([num2str(deltas(max_dir(iC,:))) '- ' sig_str])
%         end
    end
end


%% summary cell figures
figure;
for idel = 1:nDelta
    ind_cells = intersect(OSI_ind,find(max_dir(:,1) == idel));
    for ioff = 1:noff_all
        if length(ind_con)== 0 || ioff<noff_all 
            ind = setdiff(find(tFramesOff == off_all(ioff)),ind_con);
        else
            ind = ind_con;
        end
        for idelta = 1:nDelta
            ind2 = intersect(ind,find(targetDelta == deltas(idelta)));
            subplot(noff_all,nDelta,((ioff-1)*nDelta)+idel)
            errorbar(deltas(idelta),mean(nanmean(mean(data_dfof(resp_win,ind_cells,2,ind2),1)-mean(data_dfof(base_win,ind_cells,2,ind2),1),4),2),nanstd(nanmean(mean(data_dfof(resp_win,ind_cells,2,ind2),1)-mean(data_dfof(base_win,ind_cells,2,ind2),1),4),[],2)./sqrt(length(ind_cells)),'ok')
            hold on
        end
        if ioff == 1
            title(['Delta ' num2str(deltas(idel)) ' deg- ' num2str(length(ind_cells)) ' cells'])
        end
        ylim([-0.05 0.3])
    end
end
suptitle([mouse ' ' date]) 
print(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_testResp_byDelta_byInt.pdf']),'-dpdf')

figure;
g = grays(noff_all);
[n n2] = subplotn(nDelta);
for idel = 1:nDelta
    ind_cells = find(max_dir == idel);
    for ioff = 1:noff_all
        if length(ind_con)== 0 || ioff<noff_all 
            ind = setdiff(find(tFramesOff == off_all(ioff)),ind_con);
        else
            ind = ind_con;
        end
        for idelta = 1:nDelta
            ind2 = intersect(ind,find(targetDelta == deltas(idelta)));
            subplot(n,n2,idel)
            errorbar(deltas(idelta),mean(nanmean(mean(data_dfof(resp_win,ind_cells,2,ind2),1)-mean(data_dfof(base_win,ind_cells,2,ind2),1),4),2),nanstd(nanmean(mean(data_dfof(resp_win,ind_cells,2,ind2),1)-mean(data_dfof(base_win,ind_cells,2,ind2),1),4),[],2)./sqrt(length(ind_cells)),'o', 'Color', g(ioff,:))
            hold on
        end
        if ioff == 1
            title(['Delta ' num2str(deltas(idel)) ' deg- ' num2str(length(ind_cells)) ' cells'])
        end
        ylim([-0.05 0.3])
    end
end
suptitle([mouse ' ' date]) 
print(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_testResp_byDelta_byInt_overlay.pdf']),'-dpdf', '-fillpage')

figure;
g = grays(noff_all);
[n n2] = subplotn(nDelta);
for idelta = 1:nDelta
    ind_del = find(targetDelta == deltas(idelta));
    for ioff = 1:noff_all
        if length(ind_con)== 0 || ioff<noff_all 
            ind = setdiff(find(tFramesOff == off_all(ioff)),ind_con);
        else
            ind = ind_con;
        end
        for idel = 1:nDelta
            ind_cells = find(max_dir == idel);
            ind2 = intersect(ind_del,ind);
            subplot(n,n2,idelta)
            errorbar(deltas(idel),mean(nanmean(mean(data_dfof(resp_win,ind_cells,2,ind2),1)-mean(data_dfof(base_win,ind_cells,2,ind2),1),4),2),nanstd(nanmean(mean(data_dfof(resp_win,ind_cells,2,ind2),1)-mean(data_dfof(base_win,ind_cells,2,ind2),1),4),[],2)./sqrt(length(ind_cells)),'o', 'Color', g(ioff,:))
            hold on
        end
    end
    ylim([-0.05 0.3])
    xlabel('Cell preference group')
    title(['Delta ' num2str(deltas(idelta)) ' deg change'])
end
suptitle([mouse ' ' date]) 
print(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_testResp_byDeltaGroup_byInt_overlay.pdf']),'-dpdf','-fillpage')

data_dfof_resp = squeeze(mean(data_dfof(resp_win,:,:,:),1)-mean(data_dfof(base_win,:,:,:),1));
if length(ind_con)== 0
    ind_off = find(tFramesOff == off_all(end));
else
    ind_off = ind_con;
end
data_dfof_resp_max = zeros(nCells,1);
for idel = 1:nDelta
    ind_del = find(targetDelta == deltas(idel));
    ind_cells = find(max_dir == idel);
    ind = intersect(ind_off,ind_del);
    data_dfof_resp_max(ind_cells,1) = squeeze(nanmean(data_dfof_resp(ind_cells,2,ind),3));
end
data_dfof_resp_norm = bsxfun(@rdivide,squeeze(data_dfof_resp(:,2,:)), data_dfof_resp_max);

    
figure;
ori_resp = zeros(nDelta,noff_all,nDelta);
g = grays(noff_all);
[n n2] = subplotn(nDelta);
for idelta = 1:nDelta
    ind_del = find(targetDelta == deltas(idelta));
    for ioff = 1:noff_all
        if length(ind_con)== 0 || ioff<noff_all 
            ind = setdiff(find(tFramesOff == off_all(ioff)),ind_con);
        else
            ind = ind_con;
        end
        for idel = 1:nDelta
            ind_cells = find(max_dir == idel);
            ind2 = intersect(ind_del,ind);
            subplot(n,n2,idelta)
            errorbar(deltas(idel),mean(nanmean(data_dfof_resp_norm(ind_cells,ind2),2),1),nanstd(nanmean(data_dfof_resp_norm(ind_cells,ind2),2),[],1)./sqrt(length(ind_cells)),'o', 'Color', g(ioff,:))
            hold on
            ori_resp(idel,ioff,idelta) = mean(nanmean(data_dfof_resp_norm(ind_cells,ind2),2),1);
        end
    end
    ylim([-0.2 1.5])
    xlabel('Cell preference group')
    title(['Delta ' num2str(deltas(idelta)) ' deg change'])
end
suptitle([mouse ' ' date]) 
print(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_testResp_byDeltaGroup_byInt_overlay_norm.pdf']),'-dpdf','-fillpage')

%ori_resp = (nDelta,noff,nDelta), where the first column is the cell group,
%and the last column is the actual orientation- values are normalized so
%that max(ori_resp,[],3) = ones
ori_resp_thresh = ori_resp.*10;
ori_resp_thresh(find(ori_resp_thresh<0)) = 0.001;
loglike = zeros(nDelta,nDelta,nDelta,noff_all);
for idel = 1:nDelta
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            loglike(:,idel,idelta,ioff) = squeeze(log(ori_resp_thresh(idel,ioff,:)).*ori_resp_thresh(idel,noff_all,idelta));
        end
    end
end

figure;
x = get(groot,'DefaultAxesColorOrder');
maxloglike = zeros(nDelta,noff_all);
for idelta = 1:nDelta
    subplot(n,n2,idelta)
    for ioff = 1:noff_all
        h(ioff) = plot(deltas,nansum(loglike(:,:,idelta,ioff),2));
        hold on
        [max_val maxloglike(idelta,ioff)] = max(nansum(loglike(:,:,idelta,ioff),2),[],1);
    end
    if idelta == 1
        legend(num2str(off_all.*frameRateHz))
    end
    for ioff = 1:noff_all
        c = get(h(ioff), 'Color');
        plot(deltas(maxloglike(idelta,ioff)), 100,'x','Color', x(ioff,:));
        hold on
    end
    xlabel('Cell preference group')
    ylabel('Log Likelihood')
    title([num2str(deltas(idelta)) ' deg change'])
end
print(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_loglikelihood_byDelta_byInt.pdf']),'-dpdf','-fillpage')

figure; 
for ioff = 1:noff_all
    subplot(2,2,ioff)
    for idelta = 1:nDelta
        scatter(deltas(idelta),deltas(maxloglike(idelta,ioff)),'ok')
        hold on
    end
    axis square
    refline(1,0)
    title([num2str(off_all(ioff)*frameRateHz) ' ms ISI'])
    xlabel('Ori- actual')
    ylabel('Ori- max log likelihood')
end
print(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str 'scatter_maxloglike_byInt.pdf']),'-dpdf','-fillpage')
