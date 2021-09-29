clear all;
close all;
clc;
doRedChannel = 0;
ds = 'F0F1F2_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
fnout = fullfile(LG_base,'Analysis','2P');
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'PhaseRev');
svName = 'F0F1F2';
area = 'V1';
driver = 'SOM';
all_area = 0;
con = 1;
dirF0_all = [];
dirF1_all = [];
prF1_all = [];
prF2_all = [];
dfof_dir_avg_all = [];
max_ind_dir = [];
dfof_pr_avg_all = [];
max_ind_pr = [];

stim_DSI_all = [];
stim_OSI_all = [];
k1_all = [];

mouse_list = [];
    
for iexp = 1:nexp
    if strcmp(cell2mat(expt(iexp).img_loc), area) & strcmp(cell2mat(expt(iexp).driver), driver) & expt(iexp).dirCon == con
        mouse = expt(iexp).mouse;
        mouse_list = strvcat(mouse_list, mouse);
        date = expt(iexp).date;
        datemouse = [date '_' mouse];
        dirFolder = expt(iexp).dirFolder;
        nrun = length(dirFolder);
        dir_run_str = catRunName(cell2mat(dirFolder), nrun);
        datemousedirrun = [date '_' mouse '_' dir_run_str];
        prFolder = expt(iexp).prFolder;
        nrun = length(prFolder);
        pr_run_str = catRunName(cell2mat(prFolder), nrun);
        datemouseprrun = [date '_' mouse '_' pr_run_str];
        fprintf([mouse ' ' date '\n'])

        
        dirdata = load(fullfile(fnout, datemouse, datemousedirrun, [datemousedirrun '_dirAnalysis.mat']));
        prdata = load(fullfile(fnout, datemouse, datemouseprrun, [datemouseprrun '_f1f2.mat']));
        
        dirF0_all = [dirF0_all dirdata.f0];
        dirF1_all = [dirF1_all dirdata.f1];
        prF1_all = [prF1_all prdata.f1];
        prF2_all = [prF2_all prdata.f2];

        dirTCs = load(fullfile(fnout, datemouse, datemousedirrun, [datemousedirrun '_TCs.mat']));
        load(fullfile(fnout, datemouse, datemousedirrun, [datemousedirrun '_input.mat']));
        dir_input = input;
        nOff = dir_input.nScansOff;
        nOn = dir_input.nScansOn;
        tDir = celleqel2mat_padded(dir_input.tGratingDirectionDeg);
        dirs = unique(tDir);
        if find(dirs>=360)
            tDir(find(tDir>=360)) = tDir(find(tDir>=360))-360;
            dirs = unique(tDir);
        end
        nDir = length(dirs);
        nTrials = length(tDir);
        nCells = size(dirTCs.npSub_tc,2);
        trial_tc = nan(nOn+nOff,nCells,nTrials);
        for it = 1:nTrials-1
            trial_tc(:,:,it) = dirTCs.npSub_tc(1+nOff/2+((it-1).*(nOn+nOff)):nOff/2+(it.*(nOn+nOff)),:);
        end
        trial_f = mean(trial_tc(1:nOff/2,:,:),1);
        trial_df = bsxfun(@minus,trial_tc,trial_f);
        trial_dfof = bsxfun(@rdivide,trial_df,trial_f);
        dfof_dir_avg = nan(nOn+nOff,nCells,nDir,2);
        for iDir = 1:nDir
            ind = find(tDir==dirs(iDir));
            dfof_dir_avg(:,:,iDir,1) = nanmean(trial_dfof(:,:,ind),3);
            dfof_dir_avg(:,:,iDir,2) = nanstd(trial_dfof(:,:,ind),[],3);
        end
        dfof_dir_avg_all = cat(2, dfof_dir_avg_all,dfof_dir_avg);
        max_ind_dir = [max_ind_dir; dirdata.max_ind];

        prTCs = load(fullfile(fnout, datemouse, datemouseprrun, [datemouseprrun '_dfofData.mat']));
        dfof_pr_avg_all = cat(2,dfof_pr_avg_all,prTCs.data_dfof_phasedir);
        max_ind_pr = [max_ind_pr; prTCs.max_ind];
        
        load(fullfile(fnout, datemouse, datemousedirrun, [datemousedirrun '_oriResp.mat']));
        k1_all = [k1_all k1_ori];
        stim_DSI_all = [stim_DSI_all stim_DSI];
        stim_OSI_all = [stim_OSI_all stim_OSI];
    end
end
resp_ind_all = find(dirF0_all>0.02 & prF1_all>0.02);

save(fullfile(summaryDir,[svName '_Con' num2str(con) '_' area '_' driver '_summary.mat']), 'stim_DSI_all', 'stim_OSI_all', 'k1_all', 'dirF0_all', 'dirF1_all', 'prF1_all', 'prF2_all','resp_ind_all');

figure;
subplot(2,2,1)
scatter(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all),prF2_all(resp_ind_all)./prF1_all(resp_ind_all))
xlabel('F1/F0')
ylabel('F2/F1')
xlim([0 0.5])
ylim([0 4])
subplot(2,2,2)
scatter(dirF1_all(resp_ind_all),prF1_all(resp_ind_all))
xlabel('drifting F1')
ylabel('contrast rev F1')
xlim([0 1])
ylim([0 1])
refline(1)
subplot(2,2,3)
histogram(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all),[0:0.1:2])
xlabel('F1/F0')
xlim([0 4])
subplot(2,2,4)
histogram(prF2_all(resp_ind_all)./prF1_all(resp_ind_all),[0:0.1:2])
xlabel('F2/F1')
xlim([0 4])
sgtitle([driver ' ' area '- n = ' num2str(length(resp_ind_all)) ' cells; ' num2str(size(mouse_list,1)) ' mice'])
print(fullfile(summaryDir,[svName '_Con' num2str(con) '_' area '_' driver '_summary.pdf']),'-dpdf','-fillpage')


figure;
cdfplot(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all))
hold on
cdfplot(prF2_all(resp_ind_all)./prF1_all(resp_ind_all))
xlim([0 2])
print(fullfile(summaryDir,[svName '_Con' num2str(con) '_' area '_' driver '_histAll.pdf']),'-dpdf','-fillpage')

figure;
subplot(2,3,1)
scatter(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all), k1_all(resp_ind_all))
xlabel('F1/F0')
ylabel('Kappa')
lm = fitlm(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all), k1_all(resp_ind_all));
title(['p = ' num2str(chop(lm.Coefficients.pValue(2),2))])
subplot(2,3,2)
scatter(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all), stim_OSI_all(resp_ind_all))
xlabel('F1/F0')
ylabel('OSI')
lm = fitlm(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all), stim_OSI_all(resp_ind_all));
title(['p = ' num2str(chop(lm.Coefficients.pValue(2),2))])
subplot(2,3,3)
scatter(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all), stim_DSI_all(resp_ind_all))
xlabel('F1/F0')
ylabel('DSI')
lm = fitlm(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all), stim_DSI_all(resp_ind_all));
title(['p = ' num2str(chop(lm.Coefficients.pValue(2),2))])
subplot(2,3,4)
scatter(prF2_all(resp_ind_all)./prF1_all(resp_ind_all), k1_all(resp_ind_all))
xlabel('F2/F1')
ylabel('Kappa')
lm = fitlm(prF2_all(resp_ind_all)./prF1_all(resp_ind_all), k1_all(resp_ind_all));
title(['p = ' num2str(chop(lm.Coefficients.pValue(2),2))])
subplot(2,3,5)
scatter(prF2_all(resp_ind_all)./prF1_all(resp_ind_all), stim_OSI_all(resp_ind_all))
xlabel('F2/F1')
ylabel('OSI')
lm = fitlm(prF2_all(resp_ind_all)./prF1_all(resp_ind_all), stim_OSI_all(resp_ind_all));
title(['p = ' num2str(chop(lm.Coefficients.pValue(2),2))])
subplot(2,3,6)
scatter(prF2_all(resp_ind_all)./prF1_all(resp_ind_all), stim_DSI_all(resp_ind_all))
xlabel('F2/F1')
ylabel('DSI')
lm = fitlm(prF2_all(resp_ind_all)./prF1_all(resp_ind_all), stim_DSI_all(resp_ind_all));
title(['p = ' num2str(chop(lm.Coefficients.pValue(2),2))])
sgtitle([driver ' ' area])
print(fullfile(summaryDir,[svName '_Con' num2str(con) '_' area '_' driver '_Tuning.pdf']),'-dpdf','-fillpage')

%%
tt_dir = (-nOff/2+1:nOn+nOff/2).*(1000/frame_rate);
tt_pr = prTCs.tt;
ind_simp = intersect(intersect(resp_ind_all,find(dirF1_all./dirF0_all>0.15)),find(prF2_all./prF1_all<0.3));
if ~isempty(ind_simp)
    if length(ind_simp)>5
        ind_simp = randsample(ind_simp,5);
    end
    figure;
    start = 1;
    for i = 1:length(ind_simp)
        subplot(length(ind_simp),2,start)
        plot(tt_dir,dfof_dir_avg_all(:,ind_simp(i),max_ind_dir(ind_simp(i)),1));
        title(num2str(chop(dirF1_all(ind_simp(i))./dirF0_all(ind_simp(i)),2)))
        subplot(length(ind_simp),2,start+1)
        plot(tt_pr,squeeze(dfof_pr_avg_all(:,ind_simp(i),:,max_ind_pr(ind_simp(i)))));
        title(num2str(chop(prF2_all(ind_simp(i))./prF1_all(ind_simp(i)),2)))
        start = start+2;
    end
    sgtitle([driver ' ' area '-Simple cells'])
    print(fullfile(summaryDir,[svName '_Con' num2str(con) '_' area '_' driver '_ExSimpleCells.pdf']),'-dpdf','-fillpage')
end

ind_comp = intersect(intersect(resp_ind_all,find(dirF1_all./dirF0_all<0.1)),find(prF2_all./prF1_all>0.8));
if ~isempty(ind_comp)
    if length(ind_comp)>5
        ind_comp = randsample(ind_comp,5);
    end
    figure;
    start = 1;
    for i = 1:length(ind_comp)
        subplot(length(ind_comp),2,start)
        plot(tt_dir,dfof_dir_avg_all(:,ind_comp(i),max_ind_dir(ind_comp(i)),1));
        title(num2str(chop(dirF1_all(ind_comp(i))./dirF0_all(ind_comp(i)),2)))
        subplot(length(ind_comp),2,start+1)
        plot(tt_pr,squeeze(dfof_pr_avg_all(:,ind_comp(i),:,max_ind_pr(ind_comp(i)))));
        title(num2str(chop(prF2_all(ind_comp(i))./prF1_all(ind_comp(i)),2)))
        start = start+2;
    end
    sgtitle([driver ' ' area '-Complex cells'])
    print(fullfile(summaryDir,[svName '_Con' num2str(con) '_' area '_' driver '_ExComplexCells.pdf']),'-dpdf','-fillpage')
end

% ind_lowDir_highPr = randsample(intersect(find(dirF1_all./dirF0_all<0.1),intersect(find(prF2_all./prF1_all>0.5),resp_ind_all)),5);
% figure;
% start = 1;
% for i = 1:length(ind_lowDir_highPr)
%     subplot(length(ind_lowDir_highPr),2,start)
%     plot(tt_dir,dfof_dir_avg_all(:,ind_lowDir_highPr(i),max_ind_dir(ind_lowDir_highPr(i)),1));
%     title(num2str(chop(dirF1_all(ind_lowDir_highPr(i))./dirF0_all(ind_lowDir_highPr(i)),2)))
%     subplot(length(ind_lowDir_highPr),2,start+1)
%     plot(tt_pr,squeeze(dfof_pr_avg_all(:,ind_lowDir_highPr(i),:,max_ind_pr(ind_lowDir_highPr(i)))));
%     title(num2str(chop(prF2_all(ind_lowDir_highPr(i))./prF1_all(ind_lowDir_highPr(i)),2)))
%     start = start+2;
% end
% sgtitle('Low F1/F0, High F2/F1')
% print(fullfile(summaryDir,[svName '_Con' num2str(con) '_' driver '_ExLowF1F0HighF2F1Cells.pdf']),'-dpdf','-fillpage')
% 
% ind_lowDir_lowPr = randsample(intersect(find(dirF1_all./dirF0_all<0.1),intersect(find(prF2_all./prF1_all<0.5),resp_ind_all)),5);
% figure;
% start = 1;
% for i = 1:length(ind_lowDir_lowPr)
%     subplot(length(ind_lowDir_lowPr),2,start)
%     plot(tt_dir,dfof_dir_avg_all(:,ind_lowDir_lowPr(i),max_ind_dir(ind_lowDir_lowPr(i)),1));
%     title(num2str(chop(dirF1_all(ind_lowDir_lowPr(i))./dirF0_all(ind_lowDir_lowPr(i)),2)))
%     subplot(length(ind_lowDir_lowPr),2,start+1)
%     plot(tt_pr,squeeze(dfof_pr_avg_all(:,ind_lowDir_lowPr(i),:,max_ind_pr(ind_lowDir_lowPr(i)))));
%     title(num2str(chop(prF2_all(ind_lowDir_lowPr(i))./prF1_all(ind_lowDir_lowPr(i)),2)))
%     start = start+2;
% end
% sgtitle('Low F1/F0, Low F2/F1')
% print(fullfile(summaryDir,[svName '_Con' num2str(con) '_' driver '_ExLowF1F0LowF2F1Cells.pdf']),'-dpdf','-fillpage')
% 
% ind_highDir_highPr = randsample(intersect(find(dirF1_all./dirF0_all>0.1),intersect(find(prF2_all./prF1_all>0.5),resp_ind_all)),5);
% figure;
% start = 1;
% for i = 1:length(ind_highDir_highPr)
%     subplot(length(ind_highDir_highPr),2,start)
%     plot(tt_dir,dfof_dir_avg_all(:,ind_highDir_highPr(i),max_ind_dir(ind_highDir_highPr(i)),1));
%     title(num2str(chop(dirF1_all(ind_highDir_highPr(i))./dirF0_all(ind_highDir_highPr(i)),2)))
%     subplot(length(ind_highDir_highPr),2,start+1)
%     plot(tt_pr,squeeze(dfof_pr_avg_all(:,ind_highDir_highPr(i),:,max_ind_pr(ind_highDir_highPr(i)))));
%     title(num2str(chop(prF2_all(ind_highDir_highPr(i))./prF1_all(ind_highDir_highPr(i)),2)))
%     start = start+2;
% end
% sgtitle('High F1/F0, High F2/F1')
% print(fullfile(summaryDir,[svName '_Con' num2str(con) '_' driver '_ExHighF1F0HighF2F1Cells.pdf']),'-dpdf','-fillpage')

%%
driver_list = {'SCN', 'SLC','SOM','SLC'};
area_list = {'V1', 'V1', 'V1','LM'};
label_list = {'L4','L2/3','SOM','LM'};
con = 1;
f1f0_all = [];
f2f1_all = [];
cond = [];
f1f0 = zeros(length(driver_list),2);
f2f1 = zeros(length(driver_list),2);
for ic = 1:length(driver_list)
    area = area_list{ic};
    driver = driver_list{ic};
    load(fullfile(summaryDir,[svName '_Con' num2str(con) '_' area '_' driver '_summary.mat']));
    f1f0_all = [f1f0_all dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all)];
    f2f1_all = [f2f1_all prF2_all(resp_ind_all)./prF1_all(resp_ind_all)];
    cond = [cond ic.*ones(size(prF2_all(resp_ind_all)))];
    f1f0(ic,1) = mean(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all),2);
    f1f0(ic,2) = std(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all));
    f2f1(ic,1) = mean(prF2_all(resp_ind_all)./prF1_all(resp_ind_all),2);
    f2f1(ic,2) = std(prF2_all(resp_ind_all)./prF1_all(resp_ind_all),[],2)./sqrt(length(resp_ind_all));
end

[p_f1f0 t_f1f0 s_f1f0] = anova1(f1f0_all,cond,'off');
figure;
post_f1f0 = multcompare(s_f1f0);
set(gca,'YTickLabels',label_list,'XLabel','F1/F0')
[p_f2f1 t_f2f1 s_f2f1] = anova1(f2f1_all,cond,'off');
figure;
post_f2f1 = multcompare(s_f2f1);
set(gca,'YTickLabels',label_list,'XLabel','F2/F1')

figure;
narea = length(area_list);
subplot(3,2,1)
errorbar(1:length(driver_list),f1f0(:,1),f1f0(:,2),'o')
set(gca,'XTick',1:length(driver_list),'XTickLabels',label_list)
ylabel('F1/F0');
xlim([0 narea+1])
ylim([0 0.4])
hold on
sigstar(mat2cell(post_f1f0(:,1:2),ones([1 size(post_f1f0,1)])),post_f1f0(:,end))
subplot(3,2,2)
errorbar(1:length(driver_list),f2f1(:,1),f2f1(:,2),'o')
set(gca,'XTick',1:length(driver_list),'XTickLabels',label_list)
ylabel('F2/F1');
xlim([0 narea+1])
ylim([0 0.75])
hold on
sigstar(mat2cell(post_f2f1(:,1:2),ones([1 size(post_f1f0,1)])),post_f2f1(:,end))
print(fullfile(summaryDir,[svName '_Con' num2str(con) '_allExptSummary.pdf']),'-dpdf','-fillpage')
