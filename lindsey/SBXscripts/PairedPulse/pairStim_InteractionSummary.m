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

allData = cat(2, stimOne_comp_all(:,1), allStim_comp_all);
figure;
ind_grating = intersect(find(allData(:,4)<3),find(allData(:,5)<3));
subplot(1,3,1)
errorbar(1:3, nanmean(allData(ind_grating,1:3),1), nanstd(allData(ind_grating,1:3),[],1)./sqrt(length(ind_grating)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(ind_grating,ind(1)), allData(ind_grating,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['Grating comp- n = ' num2str(nGratingCells)])
xlim([0 4])
ylim([0 0.3])
ind_nat = intersect(find(allData(:,4)>2),find(allData(:,5)>2));
subplot(1,3,2)
errorbar(1:3, nanmean(allData(ind_nat,1:3),1), nanstd(allData(ind_nat,1:3),[],1)./sqrt(length(ind_nat)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(ind_nat,ind(1)), allData(ind_nat,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['Image comp- n = ' num2str(nImageCells)])
xlim([0 4])
ylim([0 0.3])
cross_mix = [intersect(find(allData(:,4)<3),find(allData(:,5)>2)); intersect(find(allData(:,4)>2),find(allData(:,5)<3))];
subplot(1,3,3)
errorbar(1:3, nanmean(allData(cross_mix,1:3),1), nanstd(allData(cross_mix,1:3),[],1)./sqrt(length(cross_mix)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(cross_mix,ind(1)), allData(cross_mix,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['Cross comp- n = ' num2str(nCrossCells)])
xlim([0 4])
ylim([0 0.3])
print(fullfile(LG_base,'Adaptation','pairStimInteraction','pairStim_GratingImage_AvgSummary.pdf'),'-dpdf','-bestfit')

%separate by grating response
figure;
ind_grating = intersect(find(allData(:,4)==2),find(allData(:,5)==1));
subplot(2,3,1)
errorbar(1:3, nanmean(allData(ind_grating,1:3),1), nanstd(allData(ind_grating,1:3),[],1)./sqrt(length(ind_grating)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(ind_grating,ind(1)), allData(ind_grating,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['0.08 cpd grating comp- n = ' num2str(length(ind_grating))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
cross_mix = intersect(find(allData(:,4)==3),find(allData(:,5)==1));
subplot(2,3,2)
errorbar(1:3, nanmean(allData(cross_mix,1:3),1), nanstd(allData(cross_mix,1:3),[],1)./sqrt(length(cross_mix)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(cross_mix,ind(1)), allData(cross_mix,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['0.08 cpd comp deer- n = ' num2str(length(cross_mix))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
cross_mix = intersect(find(allData(:,4)==4),find(allData(:,5)==1));
subplot(2,3,3)
errorbar(1:3, nanmean(allData(cross_mix,1:3),1), nanstd(allData(cross_mix,1:3),[],1)./sqrt(length(cross_mix)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(cross_mix,ind(1)), allData(cross_mix,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['0.08 cpd comp flower- n = ' num2str(length(cross_mix))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
ind_grating = intersect(find(allData(:,4)==1),find(allData(:,5)==2));
subplot(2,3,4)
errorbar(1:3, nanmean(allData(ind_grating,1:3),1), nanstd(allData(ind_grating,1:3),[],1)./sqrt(length(ind_grating)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(ind_grating,ind(1)), allData(ind_grating,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['0.32 cpd grating comp- n = ' num2str(length(ind_grating))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
cross_mix = intersect(find(allData(:,4)==3),find(allData(:,5)==2));
subplot(2,3,5)
errorbar(1:3, nanmean(allData(cross_mix,1:3),1), nanstd(allData(cross_mix,1:3),[],1)./sqrt(length(cross_mix)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(cross_mix,ind(1)), allData(cross_mix,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['0.32 cpd comp deer- n = ' num2str(length(cross_mix))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
cross_mix = intersect(find(allData(:,4)==4),find(allData(:,5)==2));
subplot(2,3,6)
errorbar(1:3, nanmean(allData(cross_mix,1:3),1), nanstd(allData(cross_mix,1:3),[],1)./sqrt(length(cross_mix)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(cross_mix,ind(1)), allData(cross_mix,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['0.32 cpd comp flower- n = ' num2str(length(cross_mix))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
print(fullfile(LG_base,'Adaptation','pairStimInteraction','pairStim_GratingResponseSummary.pdf'),'-dpdf','-bestfit')

%separate by image response
figure;
ind_img = intersect(find(allData(:,4)==4),find(allData(:,5)==3));
subplot(2,3,1)
errorbar(1:3, nanmean(allData(ind_img,1:3),1), nanstd(allData(ind_img,1:3),[],1)./sqrt(length(ind_img)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(ind_img,ind(1)), allData(ind_img,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['deer img comp- n = ' num2str(length(ind_img))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
cross_mix = intersect(find(allData(:,4)==1),find(allData(:,5)==3));
subplot(2,3,2)
errorbar(1:3, nanmean(allData(cross_mix,1:3),1), nanstd(allData(cross_mix,1:3),[],1)./sqrt(length(cross_mix)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(cross_mix,ind(1)), allData(cross_mix,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['deer comp 0.08 cpd- n = ' num2str(length(cross_mix))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
cross_mix = intersect(find(allData(:,4)==2),find(allData(:,5)==3));
subplot(2,3,3)
errorbar(1:3, nanmean(allData(cross_mix,1:3),1), nanstd(allData(cross_mix,1:3),[],1)./sqrt(length(cross_mix)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(cross_mix,ind(1)), allData(cross_mix,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['deer comp 0.32- n = ' num2str(length(cross_mix))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
ind_img = intersect(find(allData(:,4)==3),find(allData(:,5)==4));
subplot(2,3,4)
errorbar(1:3, nanmean(allData(ind_img,1:3),1), nanstd(allData(ind_img,1:3),[],1)./sqrt(length(ind_img)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(ind_img,ind(1)), allData(ind_img,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['flower comp deer- n = ' num2str(length(ind_img))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
cross_mix = intersect(find(allData(:,4)==1),find(allData(:,5)==4));
subplot(2,3,5)
errorbar(1:3, nanmean(allData(cross_mix,1:3),1), nanstd(allData(cross_mix,1:3),[],1)./sqrt(length(cross_mix)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(cross_mix,ind(1)), allData(cross_mix,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['flower comp 0.08 cpd- n = ' num2str(length(cross_mix))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
cross_mix = intersect(find(allData(:,4)==2),find(allData(:,5)==4));
subplot(2,3,6)
errorbar(1:3, nanmean(allData(cross_mix,1:3),1), nanstd(allData(cross_mix,1:3),[],1)./sqrt(length(cross_mix)),'-o')
comp = [{[1 2]}; {[1 3]}; {[2 3]}];
for i = 1:3
    ind = comp{i};
    [h(i) p(i)] = ttest(allData(cross_mix,ind(1)), allData(cross_mix,ind(2)));
end
hold on
p = 2.*p;
sigstar(comp(find(p<0.05)),p(find(p<0.05)));
title(['flower comp 0.32 cpd- n = ' num2str(length(cross_mix))])
xlim([0.5 3.5])
ylim([0 0.3])
axis square
print(fullfile(LG_base,'Adaptation','pairStimInteraction','pairStim_ImageResponseSummary.pdf'),'-dpdf','-bestfit')

save(fullfile(LG_base,'Adaptation','pairStimInteraction','pairStim_GratingImage_Summary.mat'),'allStim_comp_all','stimOne_comp_all','nGratingCells','nImageCells','nCrossCells','nCrossCellStimPairs','expt');