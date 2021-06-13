close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir_F5= fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'CrossOri_Figures', 'CrossOri_Figure5');

ds = ['CrossOriRandDirTwoPhase_ExptList'];
eval(ds);
area = 'V1';
ind = find([expt.SF] == 0.05);

nexp = size(expt,2);
totCells = 0;
totExp = 0;
exptInd = [];
resp_any_ind = [];
resp_ind_all = [];
resp_ind_45_all = [];
plaid_SI_all = [];
plaid_SI_45_all = [];
resp_stim_prefDir_all = [];
resp_mask_prefDir_all = [];
resp_plaid_prefDir_all = [];
resp_stim_45Dir_all = [];
resp_mask_45Dir_all = [];
resp_plaid_45Dir_all = [];
avg_resp_dir_all = [];
OSI_all = [];
k1_all = [];
h_resp_all= [];

for i = 1:length(ind)
    iexp = ind(i);
    
     if strcmp(expt(iexp).img_loc,area) & (strcmp(expt(iexp).driver,'SLC') || strcmp(expt(iexp).driver,'EMX'))
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        ImgFolder = expt(iexp).coFolder;
        time = expt(iexp).coTime;
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);

        fprintf([mouse ' ' date '\n'])

        %% load data

        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']));
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriFit.mat']))

        h_resp_all = [h_resp_all; h_resp(:,:,1)]; 
        plaid_SI_all = [plaid_SI_all plaid_SI(1,:)];
        plaid_SI_45_all = [plaid_SI_45_all plaid_SI_45(1,:)];
        OSI_all = [OSI_all stim_OSI];
        k1_all = [k1_all k1_ori];
        
        resp_stim_prefDir_all = [resp_stim_prefDir_all resp_stim_prefDir];
        resp_mask_prefDir_all = [resp_mask_prefDir_all resp_mask_prefDir];
        resp_plaid_prefDir_all = [resp_plaid_prefDir_all resp_plaid_prefDir(1,:)];
        resp_stim_45Dir_all = [resp_stim_45Dir_all resp_stim_45Dir];
        resp_mask_45Dir_all = [resp_mask_45Dir_all resp_mask_45Dir];
        resp_plaid_45Dir_all = [resp_plaid_45Dir_all resp_plaid_45Dir(1,:)];
        
        avg_resp_dir_all = [avg_resp_dir_all; avg_resp_dir];
        [max_val max_dir] = max(avg_resp_dir(:,:,1,1),[],2);
        resp_ind = [];
        resp_ind_45 = [];
        for iCell = 1:nCells
            dir_mask = max_dir(iCell)+4;
            if dir_mask>nStimDir
                dir_mask = dir_mask-nStimDir;
            end
            if sum(h_resp(iCell,[max_dir(iCell) dir_mask],1),2)
                resp_ind = [resp_ind iCell];
            end
            dir_45 = max_dir(iCell)+2;
            if dir_45>nStimDir
                dir_45 = dir_45-nStimDir;
            end
            dir_mask_45 = max_dir(iCell)+6;
            if dir_mask_45>nStimDir
                dir_mask_45 = dir_mask_45-nStimDir;
            end
            if sum(h_resp(iCell,[dir_45 dir_mask_45],1),2)
                resp_ind_45 = [resp_ind_45 iCell];
            end
        end
        
        resp_any = find(sum(h_resp(:,:,1,1),2))';
        resp_any_ind = [resp_any_ind resp_any+totCells];
        resp_ind_all = [resp_ind_all resp_ind+totCells];
        resp_ind_45_all = [resp_ind_45_all resp_ind_45+totCells];
        
        exptInd = [exptInd; iexp.*ones(nCells,1)];
        totCells = totCells+nCells;
        totExp = totExp + 1;
     end
end

%%
mean0_MI = mean(plaid_SI_all(resp_ind_all));
sem0_MI = std(plaid_SI_all(resp_ind_all))./sqrt(length(resp_ind_all));
med0_MI = median(plaid_SI_all(resp_ind_all));
mean45_MI = mean(plaid_SI_45_all(resp_ind_45_all));
sem45_MI = std(plaid_SI_45_all(resp_ind_45_all))./sqrt(length(resp_ind_45_all));
med45_MI = median(plaid_SI_45_all(resp_ind_45_all));
[p_MI h_MI] = ranksum(plaid_SI_all(resp_ind_all),plaid_SI_45_all(resp_ind_45_all),'tail','both');
[p_MI h_MI] = signrank(plaid_SI_all(resp_ind_45_all),plaid_SI_45_all(resp_ind_45_all),'tail','both');

figure;
subplot(3,2,1)
scatter(resp_stim_prefDir_all(resp_ind_all)+resp_mask_prefDir_all(resp_ind_all),resp_plaid_prefDir_all(resp_ind_all),'ok')
xlim([0 3])
ylim([0 3])
refline(1)
xlabel('Test + Mask')
ylabel('Plaid')
title(['Pref dir: n = ' num2str(length(resp_ind_all))])
subplot(3,2,2)
scatter(resp_stim_45Dir_all(resp_ind_45_all)+resp_mask_45Dir_all(resp_ind_45_all),resp_plaid_45Dir_all(resp_ind_45_all),'ok')
xlim([0 3])
ylim([0 3])
refline(1)
xlabel('Test + Mask')
ylabel('Plaid')
title(['45 deg from pref: n = ' num2str(length(resp_stim_45Dir_all(resp_ind_45_all)))])
subplot(3,2,3)
histogram(plaid_SI_all(resp_ind_all),-1:0.2:1)
xlim([-1 1])
xlabel('Masking index')
ylabel('Number of cells')
subplot(3,2,4)
histogram(plaid_SI_45_all(resp_ind_45_all),-1:0.2:1)
xlim([-1 1])
xlabel('Masking index')
ylabel('Number of cells')
subplot(3,2,5)
cdfplot(plaid_SI_all(resp_ind_all))
xlim([-1 1])
xlabel('Masking index')
ylabel('Fraction of cells')
title('')
subplot(3,2,6)
cdfplot(plaid_SI_45_all(resp_ind_45_all))
xlim([-1 1])
xlabel('Masking index')
ylabel('Fraction of cells')
title('')
expts = unique(exptInd);
nexp_area = length(expts);
mouse_list = [{expt.mouse}];
mice = unique(mouse_list(expts));
nmice = length(mice);
suptitle([num2str(nexp_area) ' expts; ' num2str(nmice) ' mice; ' num2str(length(resp_ind_all)) ' cells'])
print(fullfile(summaryDir_F5, 'Figure5_MI_prefV45_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_MI_prefV45_lowSF.fig'))


resp_stim_prefDir_all_rect = resp_stim_prefDir_all;
resp_stim_prefDir_all_rect(find(resp_stim_prefDir_all<0)) = 0;
resp_mask_prefDir_all_rect = resp_mask_prefDir_all;
resp_mask_prefDir_all_rect(find(resp_mask_prefDir_all<0)) = 0;
resp_stim_45Dir_all_rect = resp_stim_45Dir_all;
resp_stim_45Dir_all_rect(find(resp_stim_45Dir_all<0)) = 0;
resp_mask_45Dir_all_rect = resp_mask_45Dir_all;
resp_mask_45Dir_all_rect(find(resp_mask_45Dir_all<0)) = 0;

resp_prefDir_SI_all = abs((resp_stim_prefDir_all_rect - resp_mask_prefDir_all_rect)./(resp_stim_prefDir_all_rect + resp_mask_prefDir_all_rect));
resp_45Dir_SI_all = abs((resp_stim_45Dir_all_rect - resp_mask_45Dir_all_rect)./(resp_stim_45Dir_all_rect + resp_mask_45Dir_all_rect));

mean0_SI = mean(resp_prefDir_SI_all(resp_ind_all));
sem0_SI = std(resp_prefDir_SI_all(resp_ind_all))./sqrt(length(resp_ind_all));
med0_SI = median(resp_prefDir_SI_all(resp_ind_all));
mean45_SI = mean(resp_45Dir_SI_all(resp_ind_45_all));
sem45_SI = std(resp_45Dir_SI_all(resp_ind_45_all))./sqrt(length(resp_ind_45_all));
med45_SI = median(resp_45Dir_SI_all(resp_ind_45_all));
[p_SI h_SI] = ranksum(resp_prefDir_SI_all(resp_ind_all),resp_45Dir_SI_all(resp_ind_45_all),'tail','both');
[p_SI h_SI] = signrank(resp_prefDir_SI_all(resp_ind_45_all),resp_45Dir_SI_all(resp_ind_45_all),'tail','both');

figure;
subplot(2,2,1)
histogram(resp_prefDir_SI_all(resp_ind_all),0:0.1:1)
xlim([0 1])
xlabel('Masking index')
ylabel('Number of cells')
title(['Pref dir: n = ' num2str(length(resp_ind_all))])
subplot(2,2,2)
histogram(resp_45Dir_SI_all(resp_ind_45_all),0:0.1:1)
xlim([0 1])
xlabel('Masking index')
ylabel('Number of cells')
title(['45 deg from pref: n = ' num2str(length(resp_stim_45Dir_all(resp_ind_45_all)))])
subplot(2,2,3)
cdfplot(resp_prefDir_SI_all(resp_ind_all))
xlim([0 1])
xlabel('Masking index')
ylabel('Fraction of cells')
title('')
subplot(2,2,4)
cdfplot(resp_45Dir_SI_all(resp_ind_45_all))
xlim([0 1])
xlabel('Masking index')
ylabel('Fraction of cells')
title('')
expts = unique(exptInd);
nexp_area = length(expts);
mouse_list = [{expt.mouse}];
mice = unique(mouse_list(expts));
nmice = length(mice);
suptitle([num2str(nexp_area) ' expts; ' num2str(nmice) ' mice; ' num2str(length(resp_ind_all)) ' cells'])
print(fullfile(summaryDir_F5, 'Figure5_SI_prefV45_lowSF.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_SI_prefV45_lowSF.fig'))

%% heatmap
avg_resp_dir_all_rect = squeeze(avg_resp_dir_all(:,:,1,:,1));
avg_resp_dir_all_rect(find(avg_resp_dir_all_rect<0)) = 0;
MI_all_dir = zeros(totCells, nStimDir);
SI_all_dir = zeros(totCells, nStimDir);
for it = 1:nStimDir
    im = it+4;
    if im>nStimDir
        im = im-nStimDir;
    end
    MI_all_dir(:,it) = (avg_resp_dir_all_rect(:,it,2)-(avg_resp_dir_all_rect(:,it,1)+avg_resp_dir_all_rect(:,im,1)))./ ...
        (avg_resp_dir_all_rect(:,it,1)+avg_resp_dir_all_rect(:,im,1)+avg_resp_dir_all_rect(:,it,2));
    ind = find((avg_resp_dir_all_rect(:,it,1)+avg_resp_dir_all_rect(:,im,1))==0);
    MI_all_dir(ind,it) = nan;
    SI_all_dir(:,it) = (avg_resp_dir_all_rect(:,it,1)-avg_resp_dir_all_rect(:,im,1))./ ...
        (avg_resp_dir_all_rect(:,im,1)+avg_resp_dir_all_rect(:,it,1));
    ind = find((avg_resp_dir_all_rect(:,it,1)+avg_resp_dir_all_rect(:,im,1))==0);
    SI_all_dir(ind,it) = nan;
end
MI_all_dir_align = zeros(size(MI_all_dir));
SI_all_dir_align = zeros(size(SI_all_dir));
h_resp_all_align = zeros(size(h_resp_all));
[max_val max_dir_all] = max(avg_resp_dir_all(:,:,1),[],2);
for iCell = 1:totCells
    it = max_dir_all(iCell);
    MI_all_dir_align(iCell,:) = circshift(MI_all_dir(iCell,:),9-it);
    SI_all_dir_align(iCell,:) = circshift(SI_all_dir(iCell,:),9-it);
    h_resp_all_align(iCell,:) = circshift(h_resp_all(iCell,:),9-it);
end
nn = zeros(1,nStimDir);
resp_ind_all_dir = cell(1,nStimDir);
orthog = circshift(1:nStimDir,-4);
for it = 1:nStimDir
	nn(it) = sum(~isnan(MI_all_dir_align(resp_ind_all,it)));
    resp_ind_all_dir{it} = find(sum(h_resp_all_align(:,[it orthog(it)]),2));
end

%%
figure;
subplot(2,2,1)
nn = sum(~isnan(MI_all_dir_align(resp_ind_all,:)));
errorbar(stimDirs-180, nanmean(MI_all_dir_align(resp_ind_all,:),1), nanstd(MI_all_dir_align(resp_ind_all,:),[],1)./sqrt(nn));
xlabel('Stim-Pref Dir')
ylabel('Masking Index')
hold on
plot(0, nanmean(MI_all_dir_align(resp_ind_all,9),1),'or');
plot(45, nanmean(MI_all_dir_align(resp_ind_all,11),1),'or');
ax1 = gca; % current axes
set(ax1,'Xtick',stimDirs(1:2:end)-180, 'Ylim', [-0.6 0.2])
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'XTickLabels',-90:45:270, 'YTick', get(ax1, 'YTick'));
xlabel(ax2,'Mask-Pref Dir');

%% binned by OSI
t = 1:nStimDir;
nOSI = 5;
OSIxDiffMI = zeros(nOSI,nStimDir,2);
OSIxDiffSI = zeros(nOSI,nStimDir,2);
OSIxDiffSI_abs = zeros(nOSI,nStimDir,2);
indn = zeros(nStimDir,nOSI);
Diff_mat = [];
OSI_mat = [];
MIall_mat = [];
SIall_mat = [];
SIabsall_mat = []; 
for it = 1:nStimDir
    [n edges bin] = histcounts(OSI_all(resp_ind_all),[0:0.2:1]);
    for i = 1:length(n)
        ind = find(bin==i);
        indn(it,i) = length(ind);
        OSIxDiffMI(i,it,1) = nanmean(MI_all_dir_align(resp_ind_all(ind),t(it)),1);
        OSIxDiffMI(i,it,2) = nanstd(MI_all_dir_align(resp_ind_all(ind),t(it)),[],1)./sqrt(length(ind));
        MIall_mat = [MIall_mat; reshape(MI_all_dir_align(resp_ind_all(ind),t(it)), [length(ind) 1])];
   
    
        OSIxDiffSI(i,it,1) = nanmean(SI_all_dir_align(resp_ind_all(ind),t(it)),1);
        OSIxDiffSI(i,it,2) = nanstd(SI_all_dir_align(resp_ind_all(ind),t(it)),[],1)./sqrt(length(ind));
        SIall_mat = [SIall_mat; reshape(SI_all_dir_align(resp_ind_all(ind),t(it)), [length(ind) 1])];
    
        OSI_mat = [OSI_mat; i.*ones(length(ind),1)];
        Diff_mat = [Diff_mat; it.*ones(length(ind),1)];
    
        OSIxDiffSI_abs(i,it,1) = nanmean(abs(SI_all_dir_align(resp_ind_all(ind),t(it))),1);
        OSIxDiffSI_abs(i,it,2) = nanstd(abs(SI_all_dir_align(resp_ind_all(ind),t(it))),[],1)./sqrt(length(ind));
        SIabsall_mat = [SIabsall_mat; reshape(abs(SI_all_dir_align(resp_ind_all(ind),t(it))), [length(ind) 1])];
    end 
end
figure;
for i = 1:length(n)
     subplot(2,2,1)
    errorbar(stimDirs(1:5),OSIxDiffMI(i,9:13,1),OSIxDiffMI(i,9:13,2), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('MI')
    ylim([-0.75 0.25])
    xlim([-10 100])
    subplot(2,2,2)
    errorbar(stimDirs(1:5),OSIxDiffSI(i,9:13,1),OSIxDiffSI(i,9:13,2), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([-1 1])
    xlim([-10 100])
    subplot(2,2,3)
    errorbar(stimDirs(1:5),OSIxDiffSI_abs(i,9:13,1),OSIxDiffSI_abs(i,9:13,2), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([0 1])
    xlim([-10 100])   
end
suptitle('First 5 directions')
print(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_byOSI.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_byOSI.fig'))


figure;
for i = 1:length(n)
     subplot(2,2,1)
    errorbar(stimDirs,circshift(OSIxDiffMI(i,:,1),-8),circshift(OSIxDiffMI(i,:,2),-8), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('MI')
    ylim([-0.75 0.25])
    xlim([-10 360])
    subplot(2,2,2)
    errorbar(stimDirs,circshift(OSIxDiffSI(i,:,1),-8),circshift(OSIxDiffSI(i,:,2),-8), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([-1 1])
    xlim([-10 360])
    subplot(2,2,3)
    errorbar(stimDirs,circshift(OSIxDiffSI_abs(i,:,1),-8),circshift(OSIxDiffSI_abs(i,:,2),-8), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([0 1])
    xlim([-10 360])   
end
suptitle('All directions')
print(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_byOSI_allDir.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_byOSI_allDir.fig'))

[p_MI] = anovan(MIall_mat,{OSI_mat,Diff_mat});
[p_SI] = anovan(SIall_mat,{OSI_mat,Diff_mat});
[p_SIabs] = anovan(SIabsall_mat,{OSI_mat,Diff_mat});

%% binned by OSI only resp for each dir
t = 1:nStimDir;
nOSI = 5;
OSIxDiffMI = zeros(nOSI,nStimDir,2);
OSIxDiffSI = zeros(nOSI,nStimDir,2);
OSIxDiffSI_abs = zeros(nOSI,nStimDir,2);
indn = zeros(nStimDir,nOSI);
Diff_mat = [];
OSI_mat = [];
MIall_mat = [];
SIall_mat = [];
SIabsall_mat = []; 
for it = 1:nStimDir
    [n edges bin] = histcounts(OSI_all(resp_ind_all_dir{t(it)}),[0:0.2:1]);
    for i = 1:length(n)
        ind = find(bin==i);
        indn(it,i) = length(ind);
        OSIxDiffMI(i,it,1) = nanmean(MI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)),1);
        OSIxDiffMI(i,it,2) = nanstd(MI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)),[],1)./sqrt(length(ind));
        MIall_mat = [MIall_mat; reshape(MI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)), [length(ind) 1])];
   
    
        OSIxDiffSI(i,it,1) = nanmean(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)),1);
        OSIxDiffSI(i,it,2) = nanstd(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)),[],1)./sqrt(length(ind));
        SIall_mat = [SIall_mat; reshape(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)), [length(ind) 1])];
    
        OSI_mat = [OSI_mat; i.*ones(length(ind),1)];
        Diff_mat = [Diff_mat; it.*ones(length(ind),1)];
    
        OSIxDiffSI_abs(i,it,1) = nanmean(abs(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it))),1);
        OSIxDiffSI_abs(i,it,2) = nanstd(abs(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it))),[],1)./sqrt(length(ind));
        SIabsall_mat = [SIabsall_mat; reshape(abs(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it))), [length(ind) 1])];
    end 
end
figure;
for i = 1:length(n)
     subplot(2,2,1)
    errorbar(stimDirs(1:5),OSIxDiffMI(i,9:13,1),OSIxDiffMI(i,9:13,2), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('MI')
    ylim([-0.75 0.25])
    xlim([-10 100])
    subplot(2,2,2)
    errorbar(stimDirs(1:5),OSIxDiffSI(i,9:13,1),OSIxDiffSI(i,9:13,2), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([-1 1])
    xlim([-10 100])
    subplot(2,2,3)
    errorbar(stimDirs(1:5),OSIxDiffSI_abs(i,9:13,1),OSIxDiffSI_abs(i,9:13,2), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([0 1])
    xlim([-10 100])   
end
suptitle('First 5 directions')
print(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_byOSI_respOnly.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_byOSI_respOnly.fig'))


figure;
for i = 1:length(n)
     subplot(2,2,1)
    errorbar(stimDirs,circshift(OSIxDiffMI(i,:,1),-8),circshift(OSIxDiffMI(i,:,2),-8), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('MI')
    ylim([-0.75 0.25])
    xlim([-10 360])
    subplot(2,2,2)
    errorbar(stimDirs,circshift(OSIxDiffSI(i,:,1),-8),circshift(OSIxDiffSI(i,:,2),-8), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([-1 1])
    xlim([-10 360])
    subplot(2,2,3)
    errorbar(stimDirs,circshift(OSIxDiffSI_abs(i,:,1),-8),circshift(OSIxDiffSI_abs(i,:,2),-8), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([0 1])
    xlim([-10 360])   
end
suptitle('All directions')
print(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_byOSI_allDir_respOnly.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_byOSI_allDir_respOnly.fig'))

[p_MI] = anovan(MIall_mat,{OSI_mat,Diff_mat});
[p_SI] = anovan(SIall_mat,{OSI_mat,Diff_mat});
[p_SIabs] = anovan(SIabsall_mat,{OSI_mat,Diff_mat});

%% binned by DSI
DSI_all = zeros(size(OSI_all));
null = circshift(1:nStimDir,-8);
for iCell = 1:totCells
[max_val max_ind] = max(avg_resp_dir_all_rect(iCell,:,1,1));
null_ind = null(max_ind);
min_val = avg_resp_dir_all_rect(iCell,null_ind,1,1);
DSI_all(iCell) = (max_val-min_val)./(max_val+min_val);
end
    
t = 1:nStimDir;
nOSI = 5;
OSIxDiffMI = zeros(nOSI,nStimDir,2);
OSIxDiffSI = zeros(nOSI,nStimDir,2);
OSIxDiffSI_abs = zeros(nOSI,nStimDir,2);
indn = zeros(nStimDir,nOSI);
Diff_mat = [];
OSI_mat = [];
MIall_mat = [];
SIall_mat = [];
SIabsall_mat = []; 
for it = 1:nStimDir
    [n edges bin] = histcounts(OSI_all(resp_ind_all_dir{t(it)}),[0:0.2:1]);
    for i = 1:length(n)
        ind = intersect(find(DSI_all>0.7),find(bin==i));
        indn(it,i) = length(ind);
        OSIxDiffMI(i,it,1) = nanmean(MI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)),1);
        OSIxDiffMI(i,it,2) = nanstd(MI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)),[],1)./sqrt(length(ind));
        MIall_mat = [MIall_mat; reshape(MI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)), [length(ind) 1])];
   
    
        OSIxDiffSI(i,it,1) = nanmean(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)),1);
        OSIxDiffSI(i,it,2) = nanstd(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)),[],1)./sqrt(length(ind));
        SIall_mat = [SIall_mat; reshape(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it)), [length(ind) 1])];
    
        OSI_mat = [OSI_mat; i.*ones(length(ind),1)];
        Diff_mat = [Diff_mat; it.*ones(length(ind),1)];
    
        OSIxDiffSI_abs(i,it,1) = nanmean(abs(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it))),1);
        OSIxDiffSI_abs(i,it,2) = nanstd(abs(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it))),[],1)./sqrt(length(ind));
        SIabsall_mat = [SIabsall_mat; reshape(abs(SI_all_dir_align(resp_ind_all_dir{t(it)}(ind),t(it))), [length(ind) 1])];
    end 
end
figure;
for i = 1:length(n)
     subplot(2,2,1)
    errorbar(stimDirs(1:5),OSIxDiffMI(i,9:13,1),OSIxDiffMI(i,9:13,2), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('MI')
    ylim([-0.75 0.25])
    xlim([-10 100])
    subplot(2,2,2)
    errorbar(stimDirs(1:5),OSIxDiffSI(i,9:13,1),OSIxDiffSI(i,9:13,2), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([-1 1])
    xlim([-10 100])
    subplot(2,2,3)
    errorbar(stimDirs(1:5),OSIxDiffSI_abs(i,9:13,1),OSIxDiffSI_abs(i,9:13,2), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([0 1])
    xlim([-10 100])   
end
suptitle('First 5 directions')
print(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_byOSI_lowDSI.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_byOSI_lowDSI.fig'))


figure;
for i = 1:length(n)
     subplot(2,2,1)
    errorbar(stimDirs,circshift(OSIxDiffMI(i,:,1),-8),circshift(OSIxDiffMI(i,:,2),-8), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('MI')
    ylim([-0.75 0.25])
    xlim([-10 360])
    subplot(2,2,2)
    errorbar(stimDirs,circshift(OSIxDiffSI(i,:,1),-8),circshift(OSIxDiffSI(i,:,2),-8), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([-1 1])
    xlim([-10 360])
    subplot(2,2,3)
    errorbar(stimDirs,circshift(OSIxDiffSI_abs(i,:,1),-8),circshift(OSIxDiffSI_abs(i,:,2),-8), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([0 1])
    xlim([-10 360])   
end
suptitle('All directions')
print(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_allDir_byOSI_lowDSI.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_stimDistvMISI_allDir_byOSI_lowDSI.fig'))


[p_MI] = anovan(MIall_mat,{OSI_mat,Diff_mat});
[p_SI] = anovan(SIall_mat,{OSI_mat,Diff_mat});
[p_SIabs] = anovan(SIabsall_mat,{OSI_mat,Diff_mat});

%%

figure;
for i = 1:length(n)
    ind = find(bin==i);
    indn(i) = length(ind);
    OSIxDiff(i,:,1) = nanmean(SI_all_dir_align(resp_ind_all(ind),:),1);
    OSIxDiff(i,:,2) = std(SI_all_dir_align(resp_ind_all(ind),:),[],1)./sqrt(length(ind));
    subplot(2,2,1)
    errorbar(stimDirs(1:5),OSIxDiff(i,9:13,1),OSIxDiff(i,9:13,2), '-o')
    hold on
    xlabel('Stim-Pref dir')
    ylabel('SI')
    ylim([-1 1])
    xlim([-10 100])
    title('First 5 directions')
    
    OSIxDiffDir(i,:,1) = nanmean(SI_all_dir_shift_avg(resp_ind_all(ind),:),1);
    OSIxDiffDir(i,:,2) = std(SI_all_dir_shift_avg(resp_ind_all(ind),:),[],1)./sqrt(length(ind));
    subplot(2,2,2)
    errorbar(stimDirs(1:7),OSIxDiffDir(i,:,1),OSIxDiffDir(i,:,2), '-o')
    hold on
    xlabel('Stim/Mask-Pref dir')
    ylabel('SI')
    ylim([0 1])
    xlim([-10 180])
    title('All directions')
    
    OSIxDiffOri(i,:,1) = nanmean(SI_all_ori_shift_avg(resp_ind_all(ind),:),1);
    OSIxDiffOri(i,:,2) = std(SI_all_ori_shift_avg(resp_ind_all(ind),:),[],1)./sqrt(length(ind));
    subplot(2,2,3)
    errorbar(stimDirs(1:3),OSIxDiffOri(i,:,1),OSIxDiffOri(i,:,2), '-o')
    hold on
    xlabel('Stim/Mask-Pref ori')
    ylabel('SI')
    ylim([0 1])
    xlim([-10 55])
    title('All orientations')
end
subplot(2,2,2)
legh = legend(num2str(chop(edges(2:end)',2)),'location','southeast');
htitle= get(legh,'Title');
set(htitle,'String','OSI')
print(fullfile(summaryDir_F5, 'Figure5_stimDistvSI_byOSI.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_stimDistvSI_byOSI.fig'))


OSIxDirDiff = zeros(length(n),5);
indn = zeros(1,length(n));
for i = 1:length(n)
    ind = find(bin==i);
    indn(i) = length(ind);
    for ii = 1:nStimDir/2-1
        OSIxDirDiff(i,ii) = nanmean(nanmean(MI_all_dir_shift{ii}(ind,:),1),2);
    end
end
subplot(2,2,2)
imagesc(OSIxDirDiff); colormap('redblue'); colorbar; clim([-0.6 0.6])
set(gca,'XTick',1:nStimDir/2,'XtickLabel',stimDirs(1:nStimDir/2),'yTick',1:length(n),'YtickLabel',edges(2:end))
xlabel('Stim/Mask-Pref Dir')
ylabel('OSI')

OSIxOriDiff = zeros(length(n),5,2);
indn = zeros(1,length(n));
for i = 1:length(n)
    ind = find(bin==i);
    indn(i) = length(ind);
    for ii = 1:5
        OSIxOriDiff(i,ii,1) = nanmean(nanmean(MI_all_ori_shift{ii}(ind,:),1),2);
        OSIxOriDiff(i,ii,2) = nanstd(nanmean(MI_all_ori_shift{ii}(ind,:),2),[],1);
    end
end
subplot(2,2,3)
imagesc(OSIxOriDiff(:,:,1)); colormap('redblue'); colorbar; clim([-0.6 0.6])
set(gca,'XTick',1:5,'XtickLabel',stimDirs(1:5),'yTick',1:length(n),'YtickLabel',edges(2:end))
xlabel('Stim/Mask-Pref Ori')
ylabel('OSI')

subplot(2,2,4)
for i = 1:length(n)
plot(stimDirs(1:5),OSIxOriDiff(i,:,1))
hold on
end
xlabel('Stim/Mask-Pref Ori')
ylabel('MI')
movegui('center')
print(fullfile(summaryDir_F5, 'Figure5_MI_distFromPref_OSI_heatmap.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_MI_distFromPref_OSI_heatmap.fig'))

SIxOriDiff = nan(length(n),5,2);
ind_n = zeros(length(n),5);
for ii = 1:5
    [n edges bin] = histcounts(nanmean(SI_all_ori_shift{ii},2),[0:0.2:1]);
    for i = 1:length(n)
        ind = find(bin==i);
        ind_n(i,ii) = length(ind);
        if length(ind>0)
            SIxOriDiff(i,ii,1) = nanmean(nanmean(MI_all_ori_shift{ii}(ind,:),1),2);
            SIxOriDiff(i,ii,2) = nanstd(nanmean(MI_all_ori_shift{ii}(ind,:),2),[],1);
        end
    end
end
figure;
subplot(2,2,1)
imagesc(SIxOriDiff(:,:,1)); colormap('redblue'); colorbar; clim([-0.6 0.6])
set(gca,'XTick',1:5,'XtickLabel',stimDirs(1:5),'yTick',1:length(n),'YtickLabel',edges(2:end))
xlabel('Stim/Mask-Pref Ori')
ylabel('SI')
subplot(2,2,2)
for i = 1:length(n)
    plot(stimDirs(1:5),SIxOriDiff(i,:,1))
    hold on
end
xlabel('Stim/Mask-Pref Ori')
ylabel('MI')
ylim([-1 0])
movegui('center')
print(fullfile(summaryDir_F5, 'Figure5_MI_distFromPref_SI_heatmap.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_MI_distFromPref_SI_heatmap.fig'))


figure;
for i = 1:length(n)
    subplot(2,2,1)
    plot(stimDirs(1:5),OSIxOriDiff(i,:,1))
    hold on
    ylim([-1 .2])
    xlabel('Stim/Mask-Pref Ori')
    ylabel('MI')
    subplot(2,2,2)
    plot(stimDirs(1:5),OSIxOriDiff(i,:,2))
    hold on
    ylim([0 .5])
    xlabel('Stim/Mask-Pref Ori')
    ylabel('MI STD')
end
legend(num2str(edges(2:end)'),'location','southeast')
movegui('center')
print(fullfile(summaryDir_F5, 'Figure5_MI_distFromPref_OSI_STD.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_MI_distFromPref_OSI_STD.fig'))



[n edges bin] = histcounts(k1_all(resp_ind_all),[0:2:10]);
k1xDiff = zeros(length(n),nStimDir);
indn = zeros(1,length(n));
for i = 1:length(n)
    ind = find(bin==i);
    indn(i) = length(ind);
    k1xDiff(i,:) = nanmean(MI_all_dir_align(resp_ind_all(ind),:),1);
end
figure; 
subplot(2,2,1)
imagesc(k1xDiff); colormap('redblue'); colorbar; clim([-0.6 0.6])
set(gca,'XTick',1:2:nStimDir,'XtickLabel',stimDirs(1:2:end)-180,'yTick',1:length(n),'YtickLabel',edges(2:end))
xlabel('Stim-Pref dir')
ylabel('k')

k1xDirDiff = zeros(length(n),5);
indn = zeros(1,length(n));
for i = 1:length(n)
    ind = find(bin==i);
    indn(i) = length(ind);
    for ii = 1:nStimDir/2-1
        k1xDirDiff(i,ii) = nanmean(nanmean(MI_all_dir_shift{ii}(ind,:),2),1);
    end
end
subplot(2,2,2)
imagesc(k1xDirDiff); colormap('redblue'); colorbar; clim([-0.6 0.6])
set(gca,'XTick',1:nStimDir/2,'XtickLabel',stimDirs(1:nStimDir/2),'yTick',1:length(n),'YtickLabel',edges(2:end))
xlabel('Stim/Mask-Pref Dir')
ylabel('k')

k1xOriDiff = zeros(length(n),5);
indn = zeros(1,length(n));
for i = 1:length(n)
    ind = find(bin==i);
    indn(i) = length(ind);
    for ii = 1:5
        k1xOriDiff(i,ii) = nanmean(nanmean(MI_all_ori_shift{ii}(ind,:),2),1);
    end
end
subplot(2,2,3)
imagesc(k1xOriDiff); colormap('redblue'); colorbar; clim([-0.6 0.6])
set(gca,'XTick',1:5,'XtickLabel',stimDirs(1:5),'yTick',1:length(n),'YtickLabel',edges(2:end))
xlabel('Stim/Mask-Pref Ori')
ylabel('k')

subplot(2,2,4)
for i = 1:length(n)
plot(stimDirs(1:5),k1xOriDiff(i,:))
hold on
end
xlabel('Stim/Mask-Pref Ori')
ylabel('MI')
movegui('center')
print(fullfile(summaryDir_F5, 'Figure5_MI_distFromPref_k1_heatmap.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_MI_distFromPref_k1_heatmap.fig'))

figure; movegui('center')
subplot(2,1,1)
scatter(OSI_all(resp_ind_all),k1_all(resp_ind_all))
xlabel('OSI')
ylabel('k')
subplot(2,1,2)
scatter(OSI_all(resp_ind_all),k1_all(resp_ind_all))
xlabel('OSI')
ylabel('k')
ylim([0 8])
print(fullfile(summaryDir_F5, 'Figure5_OSI_vs_k1.pdf'),'-dpdf','-bestfit')
savefig(fullfile(summaryDir_F5, 'Figure5_OSI_vs_k1.fig'))

