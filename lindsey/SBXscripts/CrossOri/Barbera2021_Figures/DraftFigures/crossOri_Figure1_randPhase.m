close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir_F1 = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure1');

ds = ['CrossOriRandPhase_lowSF_ExptList'];
eval(ds);

%load all experiments
nexp = size(expt,2);
expt_ind = [];
resp_ind_all = [];
respTorM_ind_all = [];
resp_cell_all = cell(3,3,4);
test_resp_all = [];
mask_resp_all = [];
plaid_resp_all = [];
plaidSI_all =[];
h_plaidSI_all = [];
resp_tc_all = [];
resp_avg_all = [];
totCells = 0;
cumCells= [];
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    fprintf([mouse ' ' date '\n'])
    ImgFolder = expt(iexp).coFolder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']));
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']));
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']));
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']),'npSub_tc');
    
    resp_ind_all = [resp_ind_all; resp_ind+totCells];
    respTorM_ind_all = [respTorM_ind_all; unique([resptest_ind; respmask_ind])+totCells];
    
    test_resp = mean(resp_cell{1,3,1},2);
    test_resp(find(test_resp<0)) = 0;
    mask_resp = mean(resp_cell{3,1,1},2);
    mask_resp(find(mask_resp<0)) = 0;
    plaid_resp = mean(resp_cell{3,3,1},2);
    plaid_resp(find(plaid_resp<0)) = 0;
    
    nCells = size(test_resp,1);
    expt_ind = [expt_ind; iexp.*ones(nCells,1)];
    totCells = totCells+nCells;
    cumCells = [cumCells totCells];
    
    test_resp_all = [test_resp_all; test_resp];
    mask_resp_all = [mask_resp_all; mask_resp];
    plaid_resp_all = [plaid_resp_all; plaid_resp];
    
    plaidSI_all = [plaidSI_all; (plaid_resp-(test_resp+mask_resp))./(plaid_resp+(test_resp+mask_resp))];
    h_plaidSI = zeros(nCells,1);
    for iC = 1:nCells
        [h_plaidSI(iC,:) p_plaid] = ttest(resp_cell{3,3,1}(iC,:),test_resp(iC)+mask_resp(iC),'dim',2,'tail','both');
    end
    h_plaidSI_all = [h_plaidSI_all; h_plaidSI];
    
    postwin_frames = frame_rate.*6;
    nCells = size(test_resp,1);
    data_trial = nan(prewin_frames+postwin_frames,nCells,nTrials);
    for iTrial = 1:nTrials-1
        data_trial(:,:,iTrial) = npSub_tc(cStimOn(iTrial)-prewin_frames:cStimOn(iTrial)+postwin_frames-1,:);
    end
    data_f = mean(data_trial(1:prewin_frames,:,:),1);
    data_df = data_trial-data_f;
    data_dfof = data_df./data_f;
    
    resp_tc = zeros(prewin_frames+postwin_frames,nCells,nMaskCon,nStimCon,nMaskPhas,1);
    resp_avg = zeros(nCells,nMaskCon,nStimCon,nMaskPhas,2);
    if iexp == 1
        resp_tc_cell = cell(nMaskCon,nStimCon,nMaskPhas);
    end
    for im = 1:nMaskCon
        ind_m = find(maskCon_all == maskCons(im));
        for is = 1:nStimCon
            ind_s = find(stimCon_all == stimCons(is));
            if im>1 & is>1
                for ip = 1:nMaskPhas
                    ind_p = find(maskPhas_all == maskPhas(ip));
                    ind_use = intersect(ind_p,intersect(ind_m,ind_s));
                    resp_tc(:,:,im,is,ip,1) = nanmean(data_dfof(:,:,ind_use),3);
                    resp_tc(:,:,im,is,ip,2) = nanstd(data_dfof(:,:,ind_use),[],3)./sqrt(length(ind_use));
                    if iexp == 1
                        resp_tc_cell{im,is,ip} = data_dfof(:,:,ind_use);
                    end
                    resp_avg(:,im,is,ip,1) = squeeze(nanmean(resp_cell{im,is,ip},2));
                    resp_avg(:,im,is,ip,2) = squeeze(nanstd(resp_cell{im,is,ip},[],2))./sqrt(size(resp_cell{im,is,ip},2));
                end
            else
                ind_use = intersect(ind_m,ind_s);
                resp_tc(:,:,im,is,1,1) = nanmean(data_dfof(:,:,ind_use),3);
                resp_tc(:,:,im,is,1,2) = nanstd(data_dfof(:,:,ind_use),[],3)./sqrt(length(ind_use));
                if iexp == 1
                    resp_tc_cell{im,is,1} = data_dfof(:,:,ind_use);
                end
                resp_avg(:,im,is,1,1) = squeeze(nanmean(resp_cell{im,is,1},2));
                resp_avg(:,im,is,1,2) = squeeze(nanstd(resp_cell{im,is,1},[],2))./sqrt(size(resp_cell{im,is,1},2));
            end
        end
    end
    
    resp_tc_all = cat(2,resp_tc_all,resp_tc);
    resp_avg_all = cat(1,resp_avg_all, resp_avg);
end

% plot scatters
ind_l = intersect(respTorM_ind_all, intersect(find(plaidSI_all<0),find(h_plaidSI_all)));
ind_h = intersect(respTorM_ind_all, intersect(find(plaidSI_all>0),find(h_plaidSI_all(:,1))));
ind_n = intersect(respTorM_ind_all, find(h_plaidSI_all==0));
grating_resp_all = [test_resp_all mask_resp_all];
max_resp_all = max(grating_resp_all,[],2);
figure;
subplot(2,2,1)
scatter(max_resp_all(respTorM_ind_all),plaid_resp_all(respTorM_ind_all),'ok');
xlabel('Grating dF/F')
ylabel('Plaid dF/F')
xlim([0 3])
ylim([0 3])
refline(1)
axis square

subplot(2,2,2)
sum_resp_all = sum(grating_resp_all,2);
scatter(sum_resp_all(ind_l),plaid_resp_all(ind_l));
hold on
scatter(sum_resp_all(ind_h),plaid_resp_all(ind_h));
scatter(sum_resp_all(ind_n),plaid_resp_all(ind_n),'k');
xlabel('Test+Mask dF/F')
ylabel('Plaid dF/F')
xlim([0 3])
ylim([0 3])
refline(1)
axis square

subplot(2,2,3)
hist(plaidSI_all(respTorM_ind_all),[-1:0.2:1])
xlabel('Selectivity index- P-(T+M)/P+(T+M)')
ylabel('Number of cells')

mouse_cell = {expt.mouse};
mice  = unique(mouse_cell);
nmice = length(mice);
suptitle([num2str(nexp) ' expts; ' num2str(nmice) ' mice; ' num2str(length(respTorM_ind_all)) ' cells'])
print(fullfile(summaryDir_F1, 'Figure1_Scatters&Hist.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F1, 'Figure1_Scatters&Hist.fig'))

tt= (-prewin_frames:postwin_frames-1).*(1000./frame_rate);
figure; 
i = 1;
for im = [1 3]
    for it = 1:5
        subplot(2,5,i)
        plot(tt,resp_tc_cell{im,end,1}(:,ind_h(1),it))
        ylim([-0.5 2.5])
        if mod(i,5) == 1
            title(['T = ' num2str(stimCons(end)) '; M = ' num2str(maskCons(im))])
        end
        i = i+1;
    end
end
suptitle([expt(iexp).mouse ' ' expt(iexp).date ' - cell #' num2str(ind_h(1))])
print(fullfile(summaryDir_F1, 'Figure1_ExampleCell_posSI_singleTrialTraces.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F1, 'Figure1_ExampleCell_posSI_singleTrialTraces.fig'))

figure; 
i = 1;
for im = [1 3]
    for it = 1:5
        subplot(2,5,i)
        plot(tt,resp_tc_cell{im,end,1}(:,ind_l(5),it))
        ylim([-0.5 2.5])
        if mod(i,5) == 1
            title(['T = ' num2str(stimCons(end)) '; M = ' num2str(maskCons(im))])
        end
        i = i+1;
    end
end
suptitle([expt(iexp).mouse ' ' expt(iexp).date ' - cell #' num2str(ind_l(5))])
print(fullfile(summaryDir_F1, 'Figure1_ExampleCell_negSI_singleTrialTraces.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F1, 'Figure1_ExampleCell_negSI_singleTrialTraces.fig'))

i = 1;
figure;
for im = 1:nMaskCon
    for is = 1:nStimCon
        subplot(3,3,i)
        shadedErrorBar(tt,resp_tc_all(:,ind_h(1),im,is,1,1),resp_tc_all(:,ind_h(1),im,is,1,2));
        i = i+1;
        ylim([-.1 1.5])
        ylabel('dF/F')
        xlabel('Time (sec)')
    end
    iexp = expt_ind(ind_h(1));
end
suptitle([expt(iexp).mouse ' ' expt(iexp).date ' - cell #' num2str(ind_h(1))])
print(fullfile(summaryDir_F1, 'Figure1_ExampleCell_posSI.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F1, 'Figure1_ExampleCell_posSI.fig'))

figure;
i = 1;
for im = 1:nMaskCon
    for is = 1:nStimCon
        subplot(3,3,i)
        shadedErrorBar(tt,resp_tc_all(:,ind_l(5),im,is,1,1),resp_tc_all(:,ind_l(5),im,is,1,2));
        i = i+1;
        ylim([-.1 1.5])
        ylabel('dF/F')
        xlabel('Time (sec)')
    end
end
suptitle([expt(iexp).mouse ' ' expt(iexp).date ' - cell #' num2str(ind_l(5))])
print(fullfile(summaryDir_F1, 'Figure1_ExampleCell_negSI.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F1, 'Figure1_ExampleCell_negSI.fig'))

figure;
subplot(1,2,1)
for im = 1:nMaskCon
    errorbar(stimCons, squeeze(resp_avg_all(ind_h(1),im,:,1,1)), squeeze(resp_avg_all(ind_h(1),im,:,1,2)),'-o')
    hold on
end
xlim([0 0.5])
ylim([-.1 1])
xlabel('Grating contrast')
ylabel('dF/F')
legend(num2str(maskCons'))
title(['cell #' num2str(ind_h(1))])    

subplot(1,2,2)
for im = 1:nMaskCon
    errorbar(stimCons, squeeze(resp_avg_all(ind_l(5),im,:,1,1)), squeeze(resp_avg_all(ind_l(5),im,:,1,2)),'-o')
    hold on
end
xlim([0 0.5])
ylim([-.1 1])
title(['cell #' num2str(ind_l(5))])
xlabel('Grating contrast')
ylabel('dF/F')
legend(num2str(maskCons'))
suptitle([expt(iexp).mouse ' ' expt(iexp).date])
print(fullfile(summaryDir_F1, 'Figure1_ExampleCells_curve.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F1, 'Figure1_ExampleCells_curve.fig'))

[x y z] = scalebarCalib(str2num(date),'16x',[],1.7);
mouse = expt(iexp).mouse;
date = expt(iexp).date;
ImgFolder = expt(iexp).coFolder;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']));
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']));
nCells = max(mask_cell(:));
stats = regionprops(mask_cell);
centroids = reshape([stats.Centroid],[2 nCells])';
figure;
imagesc(data_avg)
truesize
caxis([0 2500])
colormap(gray)
hold on
plot(centroids(ind_h(1),1), centroids(ind_h(1),2),'or','MarkerSize',10)
plot(centroids(ind_l(5),1), centroids(ind_l(5),2),'ob','MarkerSize',10)
title([date ' ' mouse '- x = ' num2str(chop(x,3)) 'um; y = ' num2str(chop(y,3)) 'um'])
print(fullfile(summaryDir_F1, 'Figure1_ExampleCells_FOV.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F1, 'Figure1_ExampleCells_FOV.fig'))
