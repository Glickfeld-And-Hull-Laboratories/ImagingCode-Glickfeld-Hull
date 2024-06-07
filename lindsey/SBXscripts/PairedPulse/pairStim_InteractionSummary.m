clear all
eval(['pairStim_GratingNatImg_ExptList'])
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey\Analysis\2P';

allStim_comp_all = [];
stimOne_comp_all = [];
nGratingCells = 0;
nImageCells = 0;
nCrossCells = 0;
nCrossCellStimPairs = 0;

nexp = size(expt,2);
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    fprintf([mouse ' ' date '\n'])
    ImgFolder = expt(iexp).ImgFolder;
    nrun = size(ImgFolder,1);
    run_str = catRunName(ImgFolder, nrun);
    load(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respComp.mat']))
    allStim_comp_all = cat(1, allStim_comp_all, allStim_comp);
    stimOne_comp_all = cat(1, stimOne_comp_all, stimOne_comp);
    nGratingCells = nGratingCells + ind_combo_n(1,2);
    nImageCells = nImageCells + ind_combo_n(3,4);
    nCrossCells = nCrossCells + length(unique([ind_combo{1,3} ind_combo{1,4} ind_combo{2,3} ind_combo{2,4}]));
    nCrossCellStimPairs = nCrossCellStimPairs+ind_combo_n(1,3)+ind_combo_n(1,4)+ind_combo_n(2,3)+ind_combo_n(2,4);
end

figure;
ind_grating = intersect(find(allStim_comp_all(:,3)<3),find(allStim_comp_all(:,4)<3));
subplot(2,3,1)
plot([allStim_comp_all(ind_grating,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(allStim_comp_all(ind_grating,1:2),1), nanstd(allStim_comp_all(ind_grating,1:2),[],1)./sqrt(length(ind_grating)),'-or')
[h p_grating] = ttest(allStim_comp_all(ind_grating,1),allStim_comp_all(ind_grating,2));
title({'Grating comp', ['n = ' num2str(nGratingCells) '; p = ' num2str(chop(p_grating,3))]})
xlim([0 3])
ylim([0 1])
ind_nat = intersect(find(allStim_comp_all(:,3)>2),find(allStim_comp_all(:,4)>2));
subplot(2,3,2)
plot([allStim_comp_all(ind_nat,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(allStim_comp_all(ind_nat,1:2),1), nanstd(allStim_comp_all(ind_nat,1:2),[],1)./sqrt(length(ind_nat)),'-or')
[h p_image] = ttest(allStim_comp_all(ind_nat,1),allStim_comp_all(ind_nat,2));
title({'Image comp', ['n = ' num2str(nImageCells) '; p = ' num2str(chop(p_image,3))]})
xlim([0 3])
ylim([0 1])
cross_mix = [intersect(find(allStim_comp_all(:,3)<3),find(allStim_comp_all(:,4)>2)); intersect(find(allStim_comp_all(:,3)>2),find(allStim_comp_all(:,4)<3))];
subplot(2,3,3)
plot([allStim_comp_all(cross_mix,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(allStim_comp_all(cross_mix,1:2),1), nanstd(allStim_comp_all(cross_mix,1:2),[],1)./sqrt(length(cross_mix)),'-or')
[h p_mix] = ttest(allStim_comp_all(cross_mix,1),allStim_comp_all(cross_mix,2));
title({'Cross comp', ['n = ' num2str(nCrossCells) '; p = ' num2str(chop(p_mix,3))]})
xlim([0 3])
ylim([0 1])
ind_grating = intersect(find(allStim_comp_all(:,3)<3),find(allStim_comp_all(:,4)<3));
subplot(2,3,4)
plot([stimOne_comp_all(ind_grating,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(stimOne_comp_all(ind_grating,1:2),1), nanstd(stimOne_comp_all(ind_grating,1:2),[],1)./sqrt(length(ind_grating)),'-or')
[h p_grating] = ttest(stimOne_comp_all(ind_grating,1),stimOne_comp_all(ind_grating,2));
title({'Grating comp', ['n = ' num2str(nGratingCells) '; p = ' num2str(chop(p_grating,3))]})
xlim([0 3])
ylim([0 1])
ind_nat = intersect(find(stimOne_comp_all(:,3)>2),find(stimOne_comp_all(:,4)>2));
subplot(2,3,5)
plot([stimOne_comp_all(ind_nat,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(stimOne_comp_all(ind_nat,1:2),1), nanstd(stimOne_comp_all(ind_nat,1:2),[],1)./sqrt(length(ind_nat)),'-or')
[h p_image] = ttest(stimOne_comp_all(ind_nat,1),stimOne_comp_all(ind_nat,2));
title({'Image comp', ['n = ' num2str(nImageCells) '; p = ' num2str(chop(p_image,3))]})
xlim([0 3])
ylim([0 1])
cross_mix = [intersect(find(stimOne_comp_all(:,3)<3),find(stimOne_comp_all(:,4)>2)); intersect(find(stimOne_comp_all(:,3)>2),find(stimOne_comp_all(:,4)<3))];
subplot(2,3,6)
plot([stimOne_comp_all(cross_mix,1:2)'],'-ok')
hold on
errorbar(1:2, nanmean(stimOne_comp_all(cross_mix,1:2),1), nanstd(stimOne_comp_all(cross_mix,1:2),[],1)./sqrt(length(cross_mix)),'-or')
[h p_mix] = ttest(stimOne_comp_all(cross_mix,1),stimOne_comp_all(cross_mix,2));
title({'Cross comp', ['n = ' num2str(nCrossCells) '; p = ' num2str(chop(p_mix,3))]})
xlim([0 3])
ylim([0 1])

print(fullfile(LG_base,'Adaptation','pairStimInteraction','pairStim_GratingImage_Summary.pdf'),'-dpdf','-bestfit')

figure;
ind_grating = intersect(find(allStim_comp_all(:,3)<3),find(allStim_comp_all(:,4)<3));
subplot(1,3,1)
errorbar(1:2, nanmean(allStim_comp_all(ind_grating,1:2),1), nanstd(allStim_comp_all(ind_grating,1:2),[],1)./sqrt(length(ind_grating)),'-o')
hold on
errorbar(1, nanmean(stimOne_comp_all(ind_grating,1),1), nanstd(stimOne_comp_all(ind_grating,1),[],1)./sqrt(length(ind_grating)),'-o')
[h p_grating] = ttest(allStim_comp_all(ind_grating,1),allStim_comp_all(ind_grating,2));
[h p_grating_one] = ttest(stimOne_comp_all(ind_grating,1),stimOne_comp_all(ind_grating,2));
title({['Grating comp- n = ' num2str(nGratingCells)] ['p = ' num2str(chop(p_grating,3)) '; p = ' num2str(chop(p_grating_one,3))]})
xlim([0 3])
ylim([0 0.3])
ind_nat = intersect(find(allStim_comp_all(:,3)>2),find(allStim_comp_all(:,4)>2));
subplot(1,3,2)
errorbar(1:2, nanmean(allStim_comp_all(ind_nat,1:2),1), nanstd(allStim_comp_all(ind_nat,1:2),[],1)./sqrt(length(ind_nat)),'-o')
hold on
errorbar(1, nanmean(stimOne_comp_all(ind_nat,1),1), nanstd(stimOne_comp_all(ind_nat,1),[],1)./sqrt(length(ind_nat)),'-o')
[h p_image] = ttest(allStim_comp_all(ind_nat,1),allStim_comp_all(ind_nat,2));
[h p_image_one] = ttest(stimOne_comp_all(ind_nat,1),stimOne_comp_all(ind_nat,2));
title({['Image comp- n = ' num2str(nImageCells)] ['p = ' num2str(chop(p_image,3)) '; p = ' num2str(chop(p_image_one,3))]})
xlim([0 3])
ylim([0 0.3])
cross_mix = [intersect(find(allStim_comp_all(:,3)<3),find(allStim_comp_all(:,4)>2)); intersect(find(allStim_comp_all(:,3)>2),find(allStim_comp_all(:,4)<3))];
subplot(1,3,3)
errorbar(1:2, nanmean(allStim_comp_all(cross_mix,1:2),1), nanstd(allStim_comp_all(cross_mix,1:2),[],1)./sqrt(length(cross_mix)),'-o')
hold on
errorbar(1, nanmean(stimOne_comp_all(cross_mix,1),1), nanstd(stimOne_comp_all(cross_mix,1),[],1)./sqrt(length(cross_mix)),'-o')
[h p_mix] = ttest(allStim_comp_all(cross_mix,1),allStim_comp_all(cross_mix,2));
[h p_mix_one] = ttest(stimOne_comp_all(cross_mix,1),stimOne_comp_all(cross_mix,2));
title({['Cross comp- n = ' num2str(nCrossCells)] ['p = ' num2str(chop(p_mix,3)) '; p = ' num2str(chop(p_mix_one,3))]})
xlim([0 3])
ylim([0 0.3])
print(fullfile(LG_base,'Adaptation','pairStimInteraction','pairStim_GratingImage_AvgSummary.pdf'),'-dpdf','-bestfit')


save(fullfile(LG_base,'Adaptation','pairStimInteraction','pairStim_GratingImage_Summary.mat'),'allStim_comp_all','stimOne_comp_all','nGratingCells','nImageCells','nCrossCells','nCrossCellStimPairs','expt');