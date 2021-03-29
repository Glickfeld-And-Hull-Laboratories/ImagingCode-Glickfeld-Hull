close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir_F8 = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure8');

ds = ['CrossOriRandPhase_ExptList'];
eval(ds);
nexp = size(expt,2);

resp_ind_all = [];
respTorM_ind_all = [];
test_resp_all = [];
mask_resp_all = [];
plaid_resp_all = [];
avg_resp_dir = [];
totCells = 0;

for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    fprintf([mouse ' ' date '\n'])
    ImgFolder = expt(iexp).coFolder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    clear resp_ind
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']));
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']));
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']));
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']),'npSub_tc');
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']),'centroid_dist');
    
    resp_ind_all = [resp_ind_all; resp_ind+totCells];
    respTorM_ind_all = [respTorM_ind_all; unique([resptest_ind; respmask_ind])+totCells];
    
    nOn = frame_rate.*4;
    postwin_frames = frame_rate.*6;
    nCells = size(data_dfof_tc,2);
    totCells = totCells+nCells;
    
    data_trial = nan(prewin_frames+postwin_frames,nCells,nTrials);
    for iTrial = 1:nTrials-1
        data_trial(:,:,iTrial) = npSub_tc(cStimOn(iTrial)-prewin_frames:cStimOn(iTrial)+postwin_frames-1,:);
    end
    data_f = mean(data_trial(1:prewin_frames,:,:),1);
    data_df = data_trial-data_f;
    data_dfof = data_df./data_f;
    resp_win = prewin_frames+5:prewin_frames+nOn;
    resp_tc = zeros(prewin_frames+postwin_frames,nCells,nMaskCon,nStimCon,nMaskPhas,1);
    resp_avg = zeros(nCells,nMaskCon,nStimCon,nMaskPhas,2);
    resp_avg_rect = zeros(nCells,nMaskCon,nStimCon,nMaskPhas,2);
    resp_cell = cell(nMaskCon,nStimCon,nMaskPhas);
    resp_cell_allphase = cell(nMaskCon,nStimCon);
    ind_n = zeros(nMaskCon,nStimCon,nMaskPhas);
    for im = 1:nMaskCon
        ind_m = find(maskCon_all == maskCons(im));
        for is = 1:nStimCon
            ind_s = find(stimCon_all == stimCons(is));
            resp_cell_allphase{im,is} = [];
            if im>1 & is>1
                for ip = 1:nMaskPhas
                    ind_p = find(maskPhas_all == maskPhas(ip));
                    ind_use = intersect(find(centroid_dist<2), intersect(ind_p,intersect(ind_m,ind_s)));
                    ind_n(im,is,ip) = length(ind_use);
                    resp_tc(:,:,im,is,ip,1) = nanmean(data_dfof(:,:,ind_use),3);
                    resp_tc(:,:,im,is,ip,2) = nanstd(data_dfof(:,:,ind_use),[],3)./sqrt(length(ind_use));
                    resp_cell{im,is,ip} = squeeze(nanmean(data_dfof(resp_win,:,ind_use),1));
                    if length(ind_use)>1
                        resp_cell_allphase{im,is} = [resp_cell_allphase{im,is} resp_cell{im,is,ip}];
                    else
                        resp_cell_allphase{im,is} = [resp_cell_allphase{im,is} resp_cell{im,is,ip}'];
                    end
                    resp_cell{im,is,ip}(find(resp_cell{im,is,ip}<0)) = 0;
                    resp_avg(:,im,is,ip,1) = nanmean(resp_cell{im,is,ip},2)./sqrt(length(ind_use));
                    resp_avg(:,im,is,ip,2) = nanstd(resp_cell{im,is,ip},[],2)./sqrt(length(ind_use));
                end
            else
                ind_use = intersect(ind_m,ind_s);
                ind_n(im,is,1) = length(ind_use);
                resp_tc(:,:,im,is,1,1) = nanmean(data_dfof(:,:,ind_use),3);
                resp_tc(:,:,im,is,1,2) = nanstd(data_dfof(:,:,ind_use),[],3)./sqrt(length(ind_use));
                resp_cell{im,is,1} = squeeze(nanmean(data_dfof(resp_win,:,ind_use),1));
                resp_cell_allphase{im,is} = resp_cell{im,is,1};
                resp_avg(:,im,is,1,1) = squeeze(nanmean(resp_cell{im,is,1},2));
                resp_avg_rect(:,im,is,1,1) = resp_avg(:,im,is,1,1);
                resp_avg_rect(find(resp_avg_rect(:,im,is,1,1)<0),im,is,1,1) = 0;
                resp_avg(:,im,is,1,2) = squeeze(nanstd(resp_cell{im,is,1},[],2))./sqrt(size(resp_cell{im,is,1},2));
            end
        end
    end
    
    test_resp = resp_avg_rect(:,1,end,1,1);
    mask_resp = resp_avg_rect(:,end,1,1,1);
    plaid_resp = squeeze(resp_avg(:,end,end,:,1));
    plaid_resp(find(plaid_resp<0)) = 0;

    test_resp_all = [test_resp_all; test_resp];
    mask_resp_all = [mask_resp_all; mask_resp];
    plaid_resp_all = [plaid_resp_all; plaid_resp];
    
    ImgFolder = expt(iexp).dirFolder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
    avg_resp_dir = [avg_resp_dir; squeeze(mean(data_dfof_dir(resp_win,:,:),1))];
end

nStimDir = size(data_dfof_dir,3);
stimDirs = 0:22.5:22.5.*(nStimDir-1);
[max_val max_dir] = max(avg_resp_dir,[],2);
pop_resp_dir = nan(nStimDir,nMaskPhas+2,2);
pop_resp_comp = nan(nStimDir,nMaskPhas+2);
resp_dir_align = nan(nStimDir,nMaskPhas+2);
resp_comp_align = nan(nStimDir,nMaskPhas+2);
ind_n = zeros(1,nStimDir);

for i = 1:nStimDir
    ind = intersect(respTorM_ind_all,find(max_dir == i));
    ind_n(1,i) = length(ind);
    if length(ind)>0
        pop_resp_dir(i,1:nMaskPhas,1) = mean(plaid_resp_all(ind,:),1);
        pop_resp_dir(i,1:nMaskPhas,2) = std(plaid_resp_all(ind,:),[],1)./sqrt(length(ind));
        pop_resp_dir(i,nMaskPhas+1,1) = mean(test_resp_all(ind,:),1);
        pop_resp_dir(i,nMaskPhas+1,2) = std(test_resp_all(ind,:),[],1)./sqrt(length(ind));
        pop_resp_dir(i,nMaskPhas+2,1) = mean(mask_resp_all(ind,:),1);
        pop_resp_dir(i,nMaskPhas+2,2) = std(mask_resp_all(ind,:),[],1)./sqrt(length(ind));
    end
end
for i = 1:nMaskPhas
    figure(1)
    subplot(3,4,i)
    errorbar(stimDirs', pop_resp_dir(:,i,1), pop_resp_dir(:,i,2), '-o')
    title(num2str(maskPhas(i)))
    ylim([0 0.4])
    hold on
    figure(2)
    subplot(3,4,i)
    polarplot(deg2rad([stimDirs stimDirs(1)])', [pop_resp_dir(:,i,1); pop_resp_dir(1,i,1)])
    title(num2str(maskPhas(i)))
    rlim([0 0.4])
    hold on
end

str = strvcat('test','mask');
for ii = 1:2
    figure(1)
    subplot(3,4,i+ii)
    errorbar(stimDirs', pop_resp_dir(:,i+ii,1), pop_resp_dir(:,i+ii,2), '-o')
    title(str(ii,:))
    ylim([0 0.4])
    hold on
    figure(2)
    subplot(3,4,i+ii)
    polarplot(deg2rad([stimDirs stimDirs(1)])', [pop_resp_dir(:,i+ii,1); pop_resp_dir(1,i+ii,1)])
    title(str(ii,:))
    rlim([0 0.4])
end

figure(1)
subplot(3,4,i+ii+1)
plot(repmat(stimDirs', [1 nMaskPhas]), pop_resp_dir(:,1:nMaskPhas,1), '-')
ylim([0 0.4])
hold on
figure(2)
subplot(3,4,i+ii+1)
polarplot(deg2rad(repmat([stimDirs stimDirs(1)]',[1 nMaskPhas])), [pop_resp_dir(:,1:nMaskPhas,1); pop_resp_dir(1,1:nMaskPhas,1)])
rlim([0 0.4])


figure(1)
suptitle(num2str(ind_n))
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_byPhase_errorbar.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_byPhase_errorbar.fig'))

figure(2)
suptitle(num2str(ind_n))
print(fullfile(summaryDir_F8, 'Figure8_populationTuning_byPhase_polar.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F8, 'Figure8_populationTuning_byPhase_polarX.fig'))

    