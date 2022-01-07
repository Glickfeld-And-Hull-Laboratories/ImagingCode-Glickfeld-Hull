clear all; close all; clc;
doRedChannel = 0;
ds = 'CrossOriRandDirTwoPhase_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandDirSummary');
str = {'hiSF','lowSF'};



area_list = strvcat('V1');%,'LM');%,'AL','RL','PM');
narea = 1;%length(area_list);
driver = 'SCN';

for a = 2
    SFs = [0.1 0.05];
    ind = find([expt.SF] == SFs(a));
for iA = 1:narea
    area = area_list(iA,:);
    Zc_all = [];
    Zp_all = [];
    totCells = 0;
    resp_ind_all = [];
    resp_ind_dir_all = [];
    resp_ind_plaid_all = [];
    mouse_list = [];
    avg_resp_dir_all = [];
    plaid_corr_all = [];
    plaid_corr_rand_all = [];
    component_all = [];
    pattern_all = [];
    
for i = 1:length(ind)
    iexp = ind(i);
    if strcmp(expt(iexp).img_loc, area) & strcmp(expt(iexp).driver,driver)
    mouse = expt(iexp).mouse;
    mouse_list = strvcat(mouse_list, mouse);
    date = expt(iexp).date;
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    fprintf([mouse ' ' date '\n'])

    %% load data

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']));
    
    fprintf(['n = ' num2str(nCells) '\n'])
  
    Zc_all = [Zc_all Zc];
    Zp_all = [Zp_all Zp];
    component_all = [component_all; component];
    pattern_all = [pattern_all; pattern];
    
    avg_resp_dir_all = cat(1,avg_resp_dir_all,avg_resp_dir);
    plaid_corr_all = [plaid_corr_all plaid_corr];
    plaid_corr_rand_all = [plaid_corr_rand_all plaid_corr_rand];
    
    resp_ind = find(sum(sum(sum(h_resp,2),3),4));
    resp_ind_dir = find(sum(h_resp(:,:,1,1),2));
    resp_ind_plaid = find(sum(sum(h_resp(:,:,:,2),2),3));

    resp_ind_all = [resp_ind_all resp_ind'+totCells];
    resp_ind_dir_all = [resp_ind_dir_all resp_ind_dir'+totCells];
    resp_ind_plaid_all = [resp_ind_plaid_all resp_ind_plaid'+totCells];

    totCells = totCells+nCells;

    end
end
save(fullfile(summaryDir,['randDirTwoPhase_Summary_' str{a} '_' area '_' driver '.mat']),'mouse_list','Zc_all','Zp_all','resp_ind_all','resp_ind_dir_all','resp_ind_plaid_all', 'plaid_corr_all','avg_resp_dir_all','component_all','stimDirs')

%%
ZcZp_diff = Zc_all-Zp_all;
ind1 = intersect(resp_ind_all,intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28)));
ind2 = intersect(resp_ind_all,intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)));
Zp_use = intersect(resp_ind_all, unique([ind1 ind2]));
ind3 = intersect(resp_ind_all,intersect(find(Zc_all(1,:)>1.28),find(Zc_all(1,:)-Zp_all(1,:)>1.28)));
ind4 = intersect(resp_ind_all,intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)));
Zc_use = unique([ind3 ind4]);
ind5 = setdiff(resp_ind_all,[ind1 ind3]);
ind6 = setdiff(resp_ind_all,[ind2 ind4]);
figure; 
movegui('center')
subplot(2,2,1)
scatter(Zc_all(1,resp_ind_all), Zp_all(1,resp_ind_all))
hold on
scatter(Zc_all(1,ind1), Zp_all(1,ind1))
scatter(Zc_all(1,ind2), Zp_all(1,ind2))
xlabel('Zc')
ylabel('Zp')
ylim([-4 8])
xlim([-4 8])
plotZcZpBorders
title([num2str(maskPhas(1)) ' deg'])
subplot(2,2,2)
scatter(Zc_all(2,resp_ind_all), Zp_all(2,resp_ind_all))
hold on
scatter(Zc_all(2,ind2), Zp_all(2,ind2))
scatter(Zc_all(2,ind1), Zp_all(2,ind1))
xlabel('Zc')
ylabel('Zp')
ylim([-4 8])
xlim([-4 8])
plotZcZpBorders
title([num2str(maskPhas(2)) ' deg'])
subplot(2,2,3)
scatter(Zc_all(1,resp_ind_all), Zc_all(2,resp_ind_all))
hold on
scatter(Zc_all(1,ind1), Zc_all(2,ind1))
scatter(Zc_all(1,ind2), Zc_all(2,ind2))
xlabel(['Zc (' num2str(maskPhas(1)) ' deg)'])
ylabel(['Zc (' num2str(maskPhas(2)) ' deg)'])
ylim([-4 8])
xlim([-4 8])
refline(1)
subplot(2,2,4)
scatter(Zp_all(1,resp_ind_all), Zp_all(2,resp_ind_all))
hold on
scatter(Zp_all(1,ind1), Zp_all(2,ind1))
scatter(Zp_all(1,ind2), Zp_all(2,ind2))
xlabel(['Zp (' num2str(maskPhas(1)) ' deg)'])
ylabel(['Zp (' num2str(maskPhas(2)) ' deg)'])
ylim([-4 8])
xlim([-4 8])
refline(1)
suptitle(['RandDir Two Phase - ' str{a} ' cells- n = ' num2str(length(resp_ind_all))])
print(fullfile(summaryDir, ['randDirTwoPhase_ZpZcSummary_' str{a} '_' area '.pdf']),'-dpdf', '-fillpage')


%%
avg_resp_dir_all_circ = cat(2,avg_resp_dir_all, avg_resp_dir_all(:,1,:,:,:));
avg_resp_dir_all_circ(find(avg_resp_dir_all_circ<0)) = 0;
avg_resp_dir_all_shift = circshift(avg_resp_dir_all,2,2);
avg_resp_dir_all_circ_shift = cat(2,avg_resp_dir_all_shift, avg_resp_dir_all_shift(:,1,:,:,:));
avg_resp_dir_all_circ_shift(find(avg_resp_dir_all_circ_shift<0)) = 0;
component_all_shift = circshift(component_all,2,2);
component_all_circ = cat(2,component_all_shift, component_all_shift(:,1));
component_all_circ(find(component_all_circ<0)) = 0;
stimDirs_circ = [stimDirs stimDirs(1)];

[max_val max_dir] = max(avg_resp_dir_all(:,:,1,1,1),[],2);
align_resp_dir = zeros(totCells, nStimDir, 2, 2);
for i = 1:totCells
    align_resp_dir(i,:,:,:) = circshift(avg_resp_dir_all(i,:,:,:,1),9-max_dir(i),2);
end
align_resp_dir_circ = cat(2,align_resp_dir, align_resp_dir(:,1,:,:));
align_resp_dir_shift = circshift(align_resp_dir,2,2);
align_resp_dir_circ_shift = cat(2,align_resp_dir_shift, align_resp_dir_shift(:,1,:,:));
figure; movegui('center');
subplot(2,2,1)
polarplot(deg2rad(stimDirs_circ),mean(align_resp_dir_circ(Zc_use,:,1,1),1))
hold on
polarplot(deg2rad(stimDirs_circ),mean(cat(1,align_resp_dir_circ_shift(intersect(resp_ind_all,ind3),:,1,2),align_resp_dir_circ_shift(intersect(resp_ind_all,ind4),:,2,2)),1))
polarplot(deg2rad(stimDirs_circ),mean(cat(1,align_resp_dir_circ_shift(intersect(resp_ind_all,ind4),:,1,2),align_resp_dir_circ_shift(intersect(resp_ind_all,ind3),:,2,2)),1))
title(['Zc cells: n = ' num2str(length(Zc_use))])
subplot(2,2,2)
polarplot(deg2rad(stimDirs_circ),mean(align_resp_dir_circ(Zp_use,:,1,1),1))
hold on
polarplot(deg2rad(stimDirs_circ),mean(cat(1,align_resp_dir_circ_shift(intersect(resp_ind_all,ind1),:,1,2),align_resp_dir_circ_shift(intersect(resp_ind_all,ind2),:,2,2)),1))
polarplot(deg2rad(stimDirs_circ),mean(cat(1,align_resp_dir_circ_shift(intersect(resp_ind_all,ind2),:,1,2),align_resp_dir_circ_shift(intersect(resp_ind_all,ind1),:,2,2)),1))
title(['Zp cells: n = ' num2str(length(Zp_use))])
subplot(2,3,4)
cdfplot(plaid_corr_all(resp_ind_all))
hold on
cdfplot(plaid_corr_rand_all(resp_ind_all))
legend({'across phases','within phase'},'location','northwest')
xlim([-1 1])
xlabel('Correlation of 0/90 deg plaids')
title('')
subplot(2,3,5)
cdfplot(Zc_all(1,resp_ind_all));
hold on
cdfplot(Zc_all(2,resp_ind_all));
xlabel('Zc')
legend({'0','90'})
xlim([-4 8])
title('')
subplot(2,3,6)
cdfplot(Zp_all(1,resp_ind_all));
hold on
cdfplot(Zp_all(2,resp_ind_all));
xlabel('Zp')
legend({'0','90'})
xlim([-4 8])
title('')
% scatter(ZcZp_diff(1,resp_ind_all),ZcZp_diff(2,resp_ind_all));
% xlabel('Zc-Zp (0 deg)')
% ylabel('Zc-Zp (90 deg)')
% xlim([-8 8])
% ylim([-8 8])
% axis square
% refline(1)
suptitle([str{a} ' cells- n = ' num2str(length(resp_ind_all))])
print(fullfile(summaryDir, ['randDirTwoPhase_CorrHist_allCells_' str{a} '_' area '.pdf']),'-dpdf', '-fillpage')

figure; movegui('center');
subplot(3,2,1)
pie([length(ind1)+length(ind2) length(ind3)+length(ind4) length(ind5)+length(ind6)],{'Zp','Zc','Zn'})
title([num2str(length(resp_ind_all)) ' cells'])
subplot(3,2,2)
ZpZp = length(intersect(ind1,ind2));
ZcZc = length(intersect(ind3,ind4));
ZcZp = length(unique([intersect(ind1,ind4), intersect(ind2,ind3)]));
ZnZn = length(intersect(ind5,ind6));
ZpZn = length(unique([intersect(ind1,ind6), intersect(ind2,ind5)]));
ZcZn = length(unique([intersect(ind3,ind6), intersect(ind4,ind5)]));
pie([ZcZc ZpZp ZnZn ZcZp ZcZn ZpZn],{'Zc-Zc','Zp-Zp','Zn-Zn','Zc-Zp','Zc-Zn','Zp-Zn'})
subplot(3,2,3)
pie([ZpZp ZcZp ZpZn],{'Zp-Zp','Zc-Zp','Zp-Zn'})
title([num2str(length(Zp_use)) ' cells'])
subplot(3,2,4)
pie([ZcZc ZcZp ZcZn],{'Zc-Zc','Zc-Zp','Zc-Zn'})
title([num2str(length(Zc_use)) ' cells'])
pZp = (length(ind1)+length(ind2))./(length(resp_ind_all).*2);
nZp = ZpZp+ZpZn+ZcZp;
pZc = (length(ind3)+length(ind4))./(length(resp_ind_all).*2);
nZc = ZcZc+ZcZn+ZcZp;
pZn = (length(ind5)+length(ind6))./(length(resp_ind_all).*2);
null_ZpZp = pZp.*nZp;
null_ZcZc = pZc.*nZc;
null_ZcZp = pZp.*nZc;
null_ZpZc = pZc.*nZp;
null_ZpZn = pZn.*nZp;
null_ZcZn = pZn.*nZc;
subplot(3,2,5)
pie([null_ZpZp null_ZpZc null_ZpZn],{'Zp-Zp','Zc-Zp','Zp-Zn'})
subplot(3,2,6)
pie([null_ZcZc null_ZcZp null_ZcZn],{'Zc-Zc','Zc-Zp','Zc-Zn'})
suptitle([str{a} ' cells'])
print(fullfile(summaryDir, ['randDirTwoPhase_Pies_' str{a} '_' area '.pdf']),'-dpdf', '-fillpage')

min_val = align_resp_dir(:,1,1,1);
DSI = (max_val-min_val)./(max_val+min_val);
align_resp_ori = mean(reshape(align_resp_dir(:,:,1,1),[totCells nStimDir/2 2]),3);
max_val = align_resp_ori(:,1);
min_val = align_resp_ori(:,5);
OSI = (max_val-min_val)./(max_val+min_val);
Zc_diff = Zc_all(1,:)-Zc_all(2,:);
Zp_diff = Zp_all(1,:)-Zp_all(2,:);
figure; movegui('center');
movegui('center')
subplot(2,2,1)
cdfplot(Zc_diff(intersect(resp_ind_all,find(OSI<0.5))));
hold on
cdfplot(Zc_diff(intersect(resp_ind_all,find(OSI>0.5))));
xlabel('Zc0-Zc90')
xlim([-5 5])
legend({'OSI<0.5','OSI>0.5'},'location','southeast')
subplot(2,2,2)
cdfplot(Zp_diff(intersect(resp_ind_all,find(OSI<0.5))));
hold on
cdfplot(Zp_diff(intersect(resp_ind_all,find(OSI>0.5))));
xlabel('Zp0-Zp90')
xlim([-5 5])
legend({'OSI<0.5','OSI>0.5'},'location','southeast')
[n edges bin] = histcounts(OSI,[0:.2:1]);
Zc_OSI = zeros(length(n),2);
Zp_OSI = zeros(length(n),2);
OSI_avg = zeros(length(n),2);
for i = 1:length(n)
    Zc_OSI(i,1) = mean(Zc_diff(intersect(resp_ind_all, find(bin==i))));
    Zc_OSI(i,2) = std(Zc_diff(intersect(resp_ind_all, find(bin==i))))./sqrt(length(intersect(resp_ind_all, find(bin==i))));
    Zp_OSI(i,1) = mean(Zp_diff(intersect(resp_ind_all, find(bin==i))));
    Zp_OSI(i,2) = std(Zp_diff(intersect(resp_ind_all, find(bin==i))))./sqrt(length(intersect(resp_ind_all, find(bin==i))));
    OSI_avg(i,1) = mean(OSI(intersect(resp_ind_all, find(bin==i))));
    OSI_avg(i,2) = std(OSI(intersect(resp_ind_all, find(bin==i))))./sqrt(length(intersect(resp_ind_all, find(bin==i))));
end
subplot(2,2,3)
errorbar(OSI_avg(:,1),Zc_OSI(:,1),Zc_OSI(:,2),Zc_OSI(:,2),OSI_avg(:,2),OSI_avg(:,2));
xlabel('OSI')
ylabel('Zc0-Zc90')
ylim([-1 1])
subplot(2,2,4)
errorbar(OSI_avg(:,1),Zp_OSI(:,1),Zp_OSI(:,2),Zp_OSI(:,2),OSI_avg(:,2),OSI_avg(:,2));
xlabel('OSI')
ylabel('Zp0-Zp90')
ylim([-1 1])

figure; movegui('center');
subplot(2,2,1)
cdfplot(Zc_diff(intersect(resp_ind_all,find(DSI<0.5))));
hold on
cdfplot(Zc_diff(intersect(resp_ind_all,find(DSI>0.5))));
xlabel('Zc0-Zc90')
xlim([-5 5])
legend({'DSI<0.5','DSI>0.5'},'location','southeast')
subplot(2,2,2)
cdfplot(Zp_diff(intersect(resp_ind_all,find(DSI<0.5))));
hold on
cdfplot(Zp_diff(intersect(resp_ind_all,find(DSI>0.5))));
xlabel('Zp0-Zp90')
xlim([-5 5])
legend({'DSI<0.5','DSI>0.5'},'location','southeast')
[n edges bin] = histcounts(DSI,[0:.2:1]);
Zc_DSI = zeros(length(n),2);
Zp_DSI = zeros(length(n),2);
DSI_avg = zeros(length(n),2);
for i = 1:length(n)
    Zc_DSI(i,1) = mean(Zc_diff(intersect(resp_ind_all, find(bin==i))));
    Zc_DSI(i,2) = std(Zc_diff(intersect(resp_ind_all, find(bin==i))))./sqrt(length(intersect(resp_ind_all, find(bin==i))));
    Zp_DSI(i,1) = mean(Zp_diff(intersect(resp_ind_all, find(bin==i))));
    Zp_DSI(i,2) = std(Zp_diff(intersect(resp_ind_all, find(bin==i))))./sqrt(length(intersect(resp_ind_all, find(bin==i))));
    DSI_avg(i,1) = mean(DSI(intersect(resp_ind_all, find(bin==i))));
    DSI_avg(i,2) = std(DSI(intersect(resp_ind_all, find(bin==i))))./sqrt(length(intersect(resp_ind_all, find(bin==i))));
end
subplot(2,2,3)
errorbar(DSI_avg(:,1),Zc_DSI(:,1),Zc_DSI(:,2),Zc_DSI(:,2),DSI_avg(:,2),DSI_avg(:,2));
xlabel('DSI')
ylabel('Zc0-Zc90')
ylim([-1 1])
subplot(2,2,4)
errorbar(DSI_avg(:,1),Zp_DSI(:,1),Zp_DSI(:,2),Zp_DSI(:,2),DSI_avg(:,2),DSI_avg(:,2));
xlabel('DSI')
ylabel('Zp0-Zp90')
ylim([-1 1])

[x Zc_use] = sort(ZcZp_diff(1,:),'descend'); 
[x Zp_use] = sort(ZcZp_diff(1,:),'ascend'); 

figure; movegui('center')
for i = 1:16
    subplot(4,4,i)
    iC = Zc_use(i);
    r_max = max([avg_resp_dir_all_circ(iC,:,1,1,1) avg_resp_dir_all_circ(iC,:,1,2,1) avg_resp_dir_all_circ(iC,:,2,2,1) component_all_circ(iC,:)],[],2);
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ(iC,:,1,1,1))
    hold on
    polarplot(deg2rad(stimDirs_circ),component_all_circ(iC,:,1,1,1))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,1,2,1))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,2,2,1))
    rlim([0 r_max])
    title({['Zc0=' num2str(chop(Zc_all(1,iC),2)) ';Zc90=' num2str(chop(Zc_all(2,iC),2))], ['Zp0=' num2str(chop(Zp_all(1,iC),2)) '; Zp90=' num2str(chop(Zp_all(2,iC),2))]})
end
suptitle('Zc cells- Blue: pattern; Red: component; Yellow: plaid 0; Purple: plaid 90')
print(fullfile(summaryDir, ['randDirTwoPhase_ZcCells_' str{a} '_' area '_' num2str(n) '.pdf']),'-dpdf', '-fillpage')

figure; movegui('center')
for i = 1:16
    subplot(4,4,i)
    iC = Zp_use(i);
    r_max = max([avg_resp_dir_all_circ(iC,:,1,1,1) avg_resp_dir_all_circ(iC,:,1,2,1) avg_resp_dir_all_circ(iC,:,2,2,1) component_all_circ(iC,:)],[],2);
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ(iC,:,1,1,1))
    hold on
    polarplot(deg2rad(stimDirs_circ),component_all_circ(iC,:,1,1,1))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,1,2,1))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,2,2,1))
    rlim([0 r_max])
    title({['Zc0=' num2str(chop(Zc_all(1,iC),2)) ';Zc90=' num2str(chop(Zc_all(2,iC),2))], ['Zp0=' num2str(chop(Zp_all(1,iC),2)) '; Zp90=' num2str(chop(Zp_all(2,iC),2))]})
end
suptitle('Zp cells- Blue: pattern; Red: component; Yellow: plaid 0; Purple: plaid 90')
print(fullfile(summaryDir, ['randDirTwoPhase_ZpCells_' str{a} '_' area '_' num2str(n) '.pdf']),'-dpdf', '-fillpage')


ZcZp_ind = unique([intersect(ind1,ind4), intersect(ind2,ind3)]);
ZcZp_ind_Zc0 = intersect(ind3,ZcZp_ind);
if length(ZcZp_ind_Zc0)<16
    n = length(ZcZp_ind_Zc0);
else
    n = 16;
end

figure; movegui('center')
for i = 1:n
    subplot(4,4,i)
    iC = ZcZp_ind_Zc0(i);
    r_max = max([avg_resp_dir_all_circ(iC,:,1,1,1) avg_resp_dir_all_circ(iC,:,1,2,1) avg_resp_dir_all_circ(iC,:,2,2,1) component_all_circ(iC,:)],[],2);
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ(iC,:,1,1,1))
    hold on
    polarplot(deg2rad(stimDirs_circ),component_all_circ(iC,:,1,1,1))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,1,2,1))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,2,2,1))
    rlim([0 r_max])
    title({['Zc0=' num2str(chop(Zc_all(1,iC),2)) ';Zc90=' num2str(chop(Zc_all(2,iC),2))], ['Zp0=' num2str(chop(Zp_all(1,iC),2)) '; Zp90=' num2str(chop(Zp_all(2,iC),2))]})
end
suptitle('Zc to Zp cells- Blue: pattern; Red: component; Yellow: plaid 0; Purple: plaid 90')
print(fullfile(summaryDir, ['randDirTwoPhase_ZcToZpCells_' str{a} '_' area '_' num2str(n) '.pdf']),'-dpdf', '-fillpage')

ZcZp_ind_Zp0 = intersect(ind1,ZcZp_ind);
if length(ZcZp_ind_Zp0)<16
    n = length(ZcZp_ind_Zp0);
else
    n = 16;
end
figure; movegui('center')
for i = 1:n
    subplot(4,4,i)
    iC = ZcZp_ind_Zp0(i);
    r_max = max([avg_resp_dir_all_circ(iC,:,1,1,1) avg_resp_dir_all_circ(iC,:,1,2,1) avg_resp_dir_all_circ(iC,:,2,2,1) component_all_circ(iC,:)],[],2);
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ(iC,:,1,1,1))
    hold on
    polarplot(deg2rad(stimDirs_circ),component_all_circ(iC,:,1,1,1))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,1,2,1))
    polarplot(deg2rad(stimDirs_circ),avg_resp_dir_all_circ_shift(iC,:,2,2,1))
    rlim([0 r_max])
    title({['Zc0=' num2str(chop(Zc_all(1,iC),2)) ';Zc90=' num2str(chop(Zc_all(2,iC),2))], ['Zp0=' num2str(chop(Zp_all(1,iC),2)) '; Zp90=' num2str(chop(Zp_all(2,iC),2))]})
end
suptitle('Zp to Zc cells- Blue: pattern; Red: component; Yellow: plaid 0; Purple: plaid 90')
print(fullfile(summaryDir, ['randDirTwoPhase_ZpToZcCells_' str{a} '_' area '_' num2str(n) '.pdf']),'-dpdf', '-fillpage')

%% population tuning
%close all
pop_resp_dir = nan(nStimDir,nStimDir,nMaskPhas,2,2);
pop_resp_comp = nan(nStimDir,nStimDir,2);
resp_dir_align = nan(nStimDir,nStimDir,nMaskPhas,2,2);
resp_comp_align = nan(nStimDir,nStimDir,2);
ind_n = zeros(1,nStimDir);
for ii = 1:2
	for i = 1:nStimDir
        ind = intersect(resp_ind_all,find(max_dir == i));
        if length(ind)>0
            if ii == 2
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,:,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,:,ii,1);
            end

            pop_resp_dir(i,:,:,ii,1) = mean(temp_resp_dir(ind,:,:),1);
            pop_resp_dir(i,:,:,ii,2) = std(temp_resp_dir(ind,:,:),[],1);

            if ii == 1
                temp_resp_dir = circshift(avg_resp_dir_all(ind,:,ii,1) + circshift(avg_resp_dir_all(ind,:,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp(i,:,1)  = mean(temp_resp_dir,1);
                pop_resp_comp(i,:,2)  = std(temp_resp_dir,1);
                ind_n(1,i) = length(ind);
            end
        end
    end
end


for ii = 1:2
    for i = 1:nStimDir
        figure(1)
        subplot(5,4,i)
        errorbar(stimDirs', pop_resp_dir(:,i,1,ii,1), pop_resp_dir(:,i,1,ii,2)./sqrt(ind_n'), '-o')
        hold on
        figure(2)
        subplot(5,4,i)
        polarplot(deg2rad([stimDirs stimDirs(1)])', [pop_resp_dir(:,i,1,ii,1); pop_resp_dir(1,i,1,ii,1)])
        hold on
        if ii == 2
            figure(1)
            subplot(5,4,i)
            errorbar(stimDirs', pop_resp_dir(:,i,2,ii,1), pop_resp_dir(:,i,2,ii,2)./sqrt(ind_n'), '-o')
            hold on
            figure(2)
            subplot(5,4,i)
            polarplot(deg2rad([stimDirs stimDirs(1)])', [pop_resp_dir(:,i,2,ii,1); pop_resp_dir(1,i,2,ii,1)])
            hold on
        end
        resp_dir_align(:,i,:,ii,:) = circshift(pop_resp_dir(:,i,:,ii,:),9-i,1);
        if ii == 1
            figure(1)
            subplot(5,4,i)
            errorbar(stimDirs, pop_resp_comp(:,i,1), pop_resp_comp(:,i,2)./sqrt(ind_n'), '-o')
            hold on
            title([num2str(stimDirs(i))])
            figure(2)
            subplot(5,4,i)
            polarplot(deg2rad([stimDirs stimDirs(1)]), [pop_resp_comp(:,i,1); pop_resp_comp(1,i,1)])
            hold on
            title([num2str(stimDirs(i))])
            resp_comp_align(:,i,:) = circshift(pop_resp_comp(:,i,:),9-i,1);
        end
    end
end

for ii = 1:2
    figure(1)
    subplot(5,4,i+1)
    errorbar(stimDirs, nanmean(resp_dir_align(:,:,1,ii,1),2), nanmean(resp_dir_align(:,:,1,ii,2),2)./sqrt(length(resp_ind_all)), '-o')
    hold on
    figure(2)
    subplot(5,4,i+1)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align(:,:,1,ii,1),2); nanmean(resp_dir_align(1,:,1,ii,1),2)])
    hold on
    if ii == 2
        figure(1)
        subplot(5,4,i+1)
        errorbar(stimDirs, nanmean(resp_dir_align(:,:,2,ii,1),2), nanmean(resp_dir_align(:,:,2,ii,2),2)./sqrt(length(resp_ind_all)), '-o')
        hold on
        figure(2)
        subplot(5,4,i+1)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align(:,:,2,ii,1),2); nanmean(resp_dir_align(1,:,2,ii,1),2)])
        hold on
    end
    if ii == 1
        figure(1)
        subplot(5,4,i+1)
        errorbar(stimDirs, nanmean(resp_comp_align(:,:,1),2),nanmean(resp_comp_align(:,:,2),2)./sqrt(length(resp_ind_all)), '-o')
        hold on
        title(['All aligned'])
        figure(2)
        subplot(5,4,i+1)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align(:,:,1),2); nanmean(resp_comp_align(1,:,1),2)])
        hold on
        title(['All- aligned'])
    end
end

figure(1)
movegui('center')
suptitle({['Blue- pattern; Red- component; Yellow- plaid0; Purple- plaid90'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir, ['randDirTwoPhase_populationTuning_errorbar' str{a} '_' area '.pdf']),'-dpdf','-bestfit')

figure(2)
movegui('center')
suptitle({['Blue- pattern; Red- component; Yellow- plaid0; Purple- plaid90'], ['Cell #s- ' num2str(ind_n)]})
print(fullfile(summaryDir, ['randDirTwoPhase_populationTuning_polar' str{a} '_' area '.pdf']),'-dpdf','-bestfit')
    
pop_resp_dir_Zc0 = nan(nStimDir,nStimDir,nMaskPhas,2,2);
pop_resp_dir_Zc90 = nan(nStimDir,nStimDir,nMaskPhas,2,2);
pop_resp_dir_Zp0 = nan(nStimDir,nStimDir,nMaskPhas,2,2);
pop_resp_dir_Zp90 = nan(nStimDir,nStimDir,nMaskPhas,2,2);
pop_resp_dir_Zc = nan(nStimDir,nStimDir,nMaskPhas,2,2);
pop_resp_dir_Zp = nan(nStimDir,nStimDir,nMaskPhas,2,2);

pop_resp_comp_Zc = nan(nStimDir,nStimDir,2);
pop_resp_comp_Zp = nan(nStimDir,nStimDir,2);
pop_resp_comp_Zc0 = nan(nStimDir,nStimDir,2);
pop_resp_comp_Zp0 = nan(nStimDir,nStimDir,2);
pop_resp_comp_Zc90 = nan(nStimDir,nStimDir,2);
pop_resp_comp_Zp90 = nan(nStimDir,nStimDir,2);

resp_dir_align_Zc0 = nan(nStimDir,nStimDir,nMaskPhas,2,2);
resp_dir_align_Zc90 = nan(nStimDir,nStimDir,nMaskPhas,2,2);
resp_dir_align_Zp0 = nan(nStimDir,nStimDir,nMaskPhas,2,2);
resp_dir_align_Zp90 = nan(nStimDir,nStimDir,nMaskPhas,2,2);
resp_dir_align_Zc = nan(nStimDir,nStimDir,nMaskPhas,2,2);
resp_dir_align_Zp = nan(nStimDir,nStimDir,nMaskPhas,2,2);

resp_comp_align_Zc = nan(nStimDir,nStimDir,2);
resp_comp_align_Zp = nan(nStimDir,nStimDir,2);
resp_comp_align_Zc0 = nan(nStimDir,nStimDir,2);
resp_comp_align_Zp0 = nan(nStimDir,nStimDir,2);
resp_comp_align_Zc90 = nan(nStimDir,nStimDir,2);
resp_comp_align_Zp90 = nan(nStimDir,nStimDir,2);

for ii = 1:2
	for i = 1:nStimDir
        ind = intersect(resp_ind_all,find(max_dir == i));
        if length(ind)>0
            if ii == 2
                temp_resp_dir = circshift(avg_resp_dir_all(:,:,:,ii,1),nStimDir/8,2);
            else
                temp_resp_dir = avg_resp_dir_all(:,:,:,ii,1);
            end

            pop_resp_dir_Zc0(i,:,:,ii,1) = nanmean(temp_resp_dir(intersect(ind,ind3),:,:),1);
            pop_resp_dir_Zc90(i,:,:,ii,1) = nanmean(temp_resp_dir(intersect(ind,ind4),:,:),1);
            pop_resp_dir_Zp0(i,:,:,ii,1) = nanmean(temp_resp_dir(intersect(ind,ind1),:,:),1);
            pop_resp_dir_Zp90(i,:,:,ii,1) = nanmean(temp_resp_dir(intersect(ind,ind2),:,:),1);
            pop_resp_dir_Zc0(i,:,:,ii,2) = nanstd(temp_resp_dir(intersect(ind,ind3),:,:),[],1);
            pop_resp_dir_Zc90(i,:,:,ii,2) = nanstd(temp_resp_dir(intersect(ind,ind4),:,:),[],1);
            pop_resp_dir_Zp0(i,:,:,ii,2) = nanstd(temp_resp_dir(intersect(ind,ind1),:,:),[],1);
            pop_resp_dir_Zp90(i,:,:,ii,2) = nanstd(temp_resp_dir(intersect(ind,ind2),:,:),[],1);
            if ii == 2
                pop_resp_dir_Zc(i,:,1,ii,1) = nanmean([temp_resp_dir(intersect(ind,ind3),:,1); temp_resp_dir(intersect(ind,ind4),:,2)],1);
                pop_resp_dir_Zc(i,:,2,ii,1) = nanmean([temp_resp_dir(intersect(ind,ind3),:,2); temp_resp_dir(intersect(ind,ind4),:,1)],1);
                pop_resp_dir_Zp(i,:,1,ii,1) = nanmean([temp_resp_dir(intersect(ind,ind1),:,1); temp_resp_dir(intersect(ind,ind2),:,2)],1);
                pop_resp_dir_Zp(i,:,2,ii,1) = nanmean([temp_resp_dir(intersect(ind,ind1),:,2); temp_resp_dir(intersect(ind,ind2),:,1)],1);
                pop_resp_dir_Zc(i,:,1,ii,2) = nanstd([temp_resp_dir(intersect(ind,ind3),:,1); temp_resp_dir(intersect(ind,ind4),:,2)],[],1);
                pop_resp_dir_Zc(i,:,2,ii,2) = nanstd([temp_resp_dir(intersect(ind,ind3),:,2); temp_resp_dir(intersect(ind,ind4),:,1)],[],1);
                pop_resp_dir_Zp(i,:,1,ii,2) = nanstd([temp_resp_dir(intersect(ind,ind1),:,1); temp_resp_dir(intersect(ind,ind2),:,2)],[],1);
                pop_resp_dir_Zp(i,:,2,ii,2) = nanstd([temp_resp_dir(intersect(ind,ind1),:,2); temp_resp_dir(intersect(ind,ind2),:,1)],[],1);
            end
            if ii == 1
                pop_resp_dir_Zc(i,:,1,ii,1) = nanmean([temp_resp_dir(intersect(ind,ind3),:,1); temp_resp_dir(intersect(ind,ind4),:,1)],1);
                pop_resp_dir_Zp(i,:,1,ii,1) = nanmean([temp_resp_dir(intersect(ind,ind1),:,1); temp_resp_dir(intersect(ind,ind2),:,1)],1);
                pop_resp_dir_Zc(i,:,1,ii,2) = nanstd([temp_resp_dir(intersect(ind,ind3),:,1); temp_resp_dir(intersect(ind,ind4),:,1)],[],1);
                pop_resp_dir_Zp(i,:,1,ii,2) = nanstd([temp_resp_dir(intersect(ind,ind1),:,1); temp_resp_dir(intersect(ind,ind2),:,1)],[],1);
                indt = intersect(ind,ind3);
                temp_resp_dir = circshift(avg_resp_dir_all(indt,:,1,ii,1) + circshift(avg_resp_dir_all(indt,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp_Zc0(i,:,1)  = nanmean(temp_resp_dir,1);
                pop_resp_comp_Zc0(i,:,2)  = nanstd(temp_resp_dir,[],1);
                indt = intersect(ind,ind1);
                temp_resp_dir = circshift(avg_resp_dir_all(indt,:,1,ii,1) + circshift(avg_resp_dir_all(indt,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp_Zp0(i,:,1)  = nanmean(temp_resp_dir,1);
                pop_resp_comp_Zp0(i,:,2)  = nanstd(temp_resp_dir,[],1);
                indt = intersect(ind,ind4);
                temp_resp_dir = circshift(avg_resp_dir_all(indt,:,1,ii,1) + circshift(avg_resp_dir_all(indt,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp_Zc90(i,:,1)  = nanmean(temp_resp_dir,1);
                pop_resp_comp_Zc90(i,:,2)  = nanstd(temp_resp_dir,[],1);
                indt = intersect(ind,ind2);
                temp_resp_dir = circshift(avg_resp_dir_all(indt,:,1,ii,1) + circshift(avg_resp_dir_all(indt,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp_Zp90(i,:,1)  = nanmean(temp_resp_dir,1);
                pop_resp_comp_Zp90(i,:,2)  = nanstd(temp_resp_dir,[],1);
                indt = intersect(ind,unique([ind3 ind4]));
                temp_resp_dir = circshift(avg_resp_dir_all(indt,:,1,ii,1) + circshift(avg_resp_dir_all(indt,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp_Zc(i,:,1)  = nanmean(temp_resp_dir,1);
                pop_resp_comp_Zc(i,:,2)  = nanstd(temp_resp_dir,[],1);
                indt = intersect(ind,unique([ind1 ind2]));
                temp_resp_dir = circshift(avg_resp_dir_all(indt,:,1,ii,1) + circshift(avg_resp_dir_all(indt,:,1,ii,1),-nStimDir/4,2),nStimDir/8,2);
                pop_resp_comp_Zp(i,:,1)  = nanmean(temp_resp_dir,1);
                pop_resp_comp_Zp(i,:,2)  = nanstd(temp_resp_dir,[],1);
            end
        end
    end
end


for ii = 1:2
    for i = 1:nStimDir
        resp_dir_align_Zc0(:,i,:,ii,:) = circshift(pop_resp_dir_Zc0(:,i,:,ii,:),9-i,1);
        resp_dir_align_Zc90(:,i,:,ii,:) = circshift(pop_resp_dir_Zc90(:,i,:,ii,:),9-i,1);
        resp_dir_align_Zc(:,i,:,ii,:) = circshift(pop_resp_dir_Zc(:,i,:,ii,:),9-i,1);
        resp_dir_align_Zp0(:,i,:,ii,:) = circshift(pop_resp_dir_Zp0(:,i,:,ii,:),9-i,1);
        resp_dir_align_Zp90(:,i,:,ii,:) = circshift(pop_resp_dir_Zp90(:,i,:,ii,:),9-i,1);
        resp_dir_align_Zp(:,i,:,ii,:) = circshift(pop_resp_dir_Zp(:,i,:,ii,:),9-i,1);
        if ii == 1
            resp_comp_align_Zc0(:,i,:) = circshift(pop_resp_comp_Zc0(:,i,:),9-i,1);
            resp_comp_align_Zc90(:,i,:) = circshift(pop_resp_comp_Zc90(:,i,:),9-i,1);
            resp_comp_align_Zc(:,i,:) = circshift(pop_resp_comp_Zc(:,i,:),9-i,1);
            resp_comp_align_Zp0(:,i,:) = circshift(pop_resp_comp_Zp0(:,i,:),9-i,1);
            resp_comp_align_Zp90(:,i,:) = circshift(pop_resp_comp_Zp90(:,i,:),9-i,1);
            resp_comp_align_Zp(:,i,:) = circshift(pop_resp_comp_Zp(:,i,:),9-i,1);
        end
    end
end


for ii = 1:2
    figure(3)
    subplot(2,3,1)
    errorbar(stimDirs, nanmean(resp_dir_align_Zc0(:,:,1,ii,1),2), nanmean(resp_dir_align_Zc0(:,:,1,ii,2),2)./sqrt(length(Zc_use)), '-o')
    hold on
    subplot(2,3,3)
    errorbar(stimDirs, nanmean(resp_dir_align_Zc(:,:,1,ii,1),2), nanmean(resp_dir_align_Zc(:,:,1,ii,2),2)./sqrt(length(Zc_use)), '-o')
    hold on
    subplot(2,3,4)
    errorbar(stimDirs, nanmean(resp_dir_align_Zp0(:,:,1,ii,1),2), nanmean(resp_dir_align_Zp0(:,:,1,ii,2),2)./sqrt(length(Zp_use)), '-o')
    hold on
    subplot(2,3,6)
    errorbar(stimDirs, nanmean(resp_dir_align_Zp(:,:,1,ii,1),2), nanmean(resp_dir_align_Zp(:,:,1,ii,2),2)./sqrt(length(Zp_use)), '-o')
    hold on
    figure(4)
    subplot(2,3,1)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc0(:,:,1,ii,1),2); nanmean(resp_dir_align_Zc0(1,:,1,ii,1),2)])
    hold on
    subplot(2,3,3)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc(:,:,1,ii,1),2); nanmean(resp_dir_align_Zc(1,:,1,ii,1),2)]) 
    hold on
    subplot(2,3,4)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp0(:,:,1,ii,1),2); nanmean(resp_dir_align_Zp0(1,:,1,ii,1),2)])
    hold on
    subplot(2,3,6)
    polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp(:,:,1,ii,1),2); nanmean(resp_dir_align_Zp(1,:,1,ii,1),2)])
    hold on
    if ii == 2
        figure(3)
        subplot(2,3,2)
        errorbar(stimDirs, nanmean(resp_dir_align_Zc90(:,:,2,ii,1),2), nanmean(resp_dir_align_Zc90(:,:,2,ii,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        subplot(2,3,5)
        errorbar(stimDirs, nanmean(resp_dir_align_Zp90(:,:,2,ii,1),2), nanmean(resp_dir_align_Zp90(:,:,2,ii,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        subplot(2,3,1)
        errorbar(stimDirs, nanmean(resp_dir_align_Zc0(:,:,2,ii,1),2), nanmean(resp_dir_align_Zc0(:,:,2,ii,2),2)./sqrt(length(Zc_use)), '-o')
        title('Zc- 0 deg')
        hold on
        subplot(2,3,2)
        errorbar(stimDirs, nanmean(resp_dir_align_Zc90(:,:,1,ii,1),2), nanmean(resp_dir_align_Zc90(:,:,1,ii,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        title('Zc- 90 deg')
        subplot(2,3,3)
        errorbar(stimDirs, nanmean(resp_dir_align_Zc(:,:,2,ii,1),2), nanmean(resp_dir_align_Zc(:,:,2,ii,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        title('Zc- all')
        subplot(2,3,4)
        errorbar(stimDirs, nanmean(resp_dir_align_Zp0(:,:,2,ii,1),2), nanmean(resp_dir_align_Zp0(:,:,2,ii,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        title('Zp- 0 deg')
        subplot(2,3,5)
        errorbar(stimDirs, nanmean(resp_dir_align_Zp90(:,:,1,ii,1),2), nanmean(resp_dir_align_Zp90(:,:,1,ii,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        title('Zp- 90 deg')
        subplot(2,3,6)
        errorbar(stimDirs, nanmean(resp_dir_align_Zp(:,:,2,ii,1),2), nanmean(resp_dir_align_Zp(:,:,2,ii,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        title('Zp- all')
        figure(4)
        hold on
        subplot(2,3,2)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc90(:,:,2,ii,1),2); nanmean(resp_dir_align_Zc90(1,:,2,ii,1),2)]) 
        hold on
        subplot(2,3,5)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp90(:,:,2,ii,1),2); nanmean(resp_dir_align_Zp90(1,:,2,ii,1),2)])
        hold on
        subplot(2,3,1)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc0(:,:,2,ii,1),2); nanmean(resp_dir_align_Zc0(1,:,2,ii,1),2)])
        hold on
        title('Zc- 0 deg')
        subplot(2,3,2)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc90(:,:,1,ii,1),2); nanmean(resp_dir_align_Zc90(1,:,1,ii,1),2)]) 
        hold on
        title('Zc- 90 deg')
        subplot(2,3,3)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc(:,:,2,ii,1),2); nanmean(resp_dir_align_Zc(1,:,2,ii,1),2)]) 
        hold on
        title('Zc- all')
        subplot(2,3,4)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp0(:,:,2,ii,1),2); nanmean(resp_dir_align_Zp0(1,:,2,ii,1),2)])
        hold on
        title('Zp- 0 deg')
        subplot(2,3,5)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp90(:,:,1,ii,1),2); nanmean(resp_dir_align_Zp90(1,:,1,ii,1),2)])
        hold on
        title('Zp- 90 deg')
        subplot(2,3,6)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp(:,:,2,ii,1),2); nanmean(resp_dir_align_Zp(1,:,2,ii,1),2)])
        hold on
        title('Zp- all')
    end
    if ii == 1
        figure(3)
        subplot(2,3,2)
        errorbar(stimDirs, nanmean(resp_dir_align_Zc90(:,:,1,ii,1),2), nanmean(resp_dir_align_Zc90(:,:,1,ii,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        subplot(2,3,5)
        errorbar(stimDirs, nanmean(resp_dir_align_Zp90(:,:,1,ii,1),2), nanmean(resp_dir_align_Zp90(:,:,1,ii,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        subplot(2,3,1)
        errorbar(stimDirs, nanmean(resp_comp_align_Zc0(:,:,1),2),nanmean(resp_comp_align_Zc0(:,:,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        ylim([0 0.4])
        subplot(2,3,2)
        errorbar(stimDirs, nanmean(resp_comp_align_Zc90(:,:,1),2), nanmean(resp_comp_align_Zc90(:,:,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        ylim([0 0.4])
        subplot(2,3,3)
        errorbar(stimDirs, nanmean(resp_comp_align_Zc(:,:,1),2), nanmean(resp_comp_align_Zc(:,:,2),2)./sqrt(length(Zc_use)), '-o')
        hold on
        ylim([0 0.4])
        subplot(2,3,4)
        errorbar(stimDirs, nanmean(resp_comp_align_Zp0(:,:,1),2), nanmean(resp_comp_align_Zp0(:,:,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        ylim([0 0.4])
        subplot(2,3,5)
        errorbar(stimDirs, nanmean(resp_comp_align_Zp90(:,:,1),2), nanmean(resp_comp_align_Zp90(:,:,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        ylim([0 0.4])
        subplot(2,3,6)
        errorbar(stimDirs, nanmean(resp_comp_align_Zp(:,:,1),2), nanmean(resp_comp_align_Zp(:,:,2),2)./sqrt(length(Zp_use)), '-o')
        hold on
        ylim([0 0.4])
        figure(4)
        hold on
        subplot(2,3,2)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zc90(:,:,1,ii,1),2); nanmean(resp_dir_align_Zc90(1,:,1,ii,1),2)]) 
        hold on
        subplot(2,3,5)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_dir_align_Zp90(:,:,1,ii,1),2); nanmean(resp_dir_align_Zp90(1,:,1,ii,1),2)])
        hold on
        subplot(2,3,1)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zc0(:,:,1),2); nanmean(resp_comp_align_Zc0(1,:,1),2)])
        hold on
        subplot(2,3,2)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zc90(:,:,1),2); nanmean(resp_comp_align_Zc90(1,:,1),2)]) 
        hold on
        subplot(2,3,3)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zc(:,:,1),2); nanmean(resp_comp_align_Zc(1,:,1),2)]) 
        hold on
        subplot(2,3,4)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zp0(:,:,1),2); nanmean(resp_comp_align_Zp0(1,:,1),2)])
        hold on
        subplot(2,3,5)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zp90(:,:,1),2); nanmean(resp_comp_align_Zp90(1,:,1),2)])
        hold on
        subplot(2,3,6)
        polarplot(deg2rad([stimDirs stimDirs(1)]), [nanmean(resp_comp_align_Zp(:,:,1),2); nanmean(resp_comp_align_Zp(1,:,1),2)])
        hold on
    end
end

figure(3)
movegui('center')
suptitle({['Top Row: Zc Cells- ' num2str(length(Zc_use))], ['Bottom Row: Zp Cells- ' num2str(length(Zp_use))]})
print(fullfile(summaryDir, ['randDirTwoPhase_populationTuning_ZpZcCells_errorbar' str{a} '_' area '.pdf']),'-dpdf','-bestfit')

figure(4)
movegui('center')
suptitle({['Top Row: Zc Cells- ' num2str(length(Zc_use))], ['Bottom Row: Zp Cells- ' num2str(length(Zp_use))]})
print(fullfile(summaryDir, ['randDirTwoPhase_populationTuning_ZpZcCells_polar' str{a} '_' area '.pdf']),'-dpdf','-bestfit')

end
end