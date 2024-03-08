close all
clear all
clc
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'Adaptation', 'SFSummary');

area_list = strvcat('V1','LM');
nDir = 6;

figure;
for iArea = 1:length(area_list)
    area = area_list(iArea,:);
    load(fullfile(summaryDir,['adaptationDirBySF_' area '_' num2str(nDir) 'Dirs.mat']))
    subplot(2,2,1)
    errorbar(sfs,norm_sf_avg(:,1),norm_sf_avg(:,2),'-o')
    hold on
    subplot(2,2,2)
    errorbar(sfs,norm_prefsf_avg(:,1),norm_prefsf_avg(:,2),'-o')
    hold on
    subplot(2,2,3)
    errorbar(oris,norm_ori_avg(:,1),norm_ori_avg(:,2),'-o')
    hold on
    subplot(2,2,4)
    errorbar(oris,norm_prefori_avg(:,1),norm_prefori_avg(:,2),'-o')
    hold on
end
subplot(2,2,1)
title('Responsive')
xlabel('SFs')
ylabel('Norm. response')
ylim([0 2])
legend(area_list)
subplot(2,2,2)
title('Preferred')
xlabel('SFs')
ylabel('Norm. response')
ylim([0 2])
subplot(2,2,3)
xlabel('Oris')
ylabel('Norm. response')
ylim([0 2])
subplot(2,2,4)
xlabel('Oris')
ylabel('Norm. response')
ylim([0 2])


