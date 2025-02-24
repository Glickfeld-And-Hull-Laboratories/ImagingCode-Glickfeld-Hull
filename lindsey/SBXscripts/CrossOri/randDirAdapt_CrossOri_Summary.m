clc; clear all; close all;
ds = 'CrossOriSingleStimRandDirAdapt_ExptList';
eval(ds)
nexp = size(expt,2);

frame_rate = 15;

Zp_all = [];
Zc_all = [];
avg_resp_dir_shift_all = [];
maxDir_all = [];
maxVal_all = [];
DSI_all = [];
h_resp_all = [];

for iexp  = 1:nexp

    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    ImgFolder = expt(iexp).coFolder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    saveLoc = expt(iexp).saveLoc;
    
    base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\';
    
    fprintf([mouse ' ' date '\n'])
    
    load(fullfile(base, saveLoc, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(base, saveLoc, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
    load(fullfile(base, saveLoc, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(base, saveLoc, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    
    maskDiff_all = celleqel2mat_padded(input.tMaskTwoGratingDirectionDeg) - celleqel2mat_padded(input.tStimTwoGratingDirectionDeg);
    maskDiffs = unique(maskDiff_all);
    nCells = size(avg_resp_dir,1);
    
    testDirs = unique(celleqel2mat_padded(input.tStimTwoGratingDirectionDeg));
    nTestDirs = length(testDirs);
    
    int = unique(diff(testDirs));
    component = squeeze(avg_resp_dir(:,:,1,:,1)+circshift(avg_resp_dir(:,:,1,:,1),-maskDiffs./int,2));
    pattern = squeeze(circshift(avg_resp_dir(:,:,1,:,1),-(maskDiffs/2)./int,2));
    
    comp_corr = zeros(2,nCells);
    patt_corr = zeros(2,nCells);
    comp_patt_corr = zeros(2,nCells);
    
    for iCell = 1:nCells
        for i = 1:2
            comp_corr(i,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,i,1),component(iCell,:,i)));
            patt_corr(i,iCell) = triu2vec(corrcoef(avg_resp_dir(iCell,:,2,i,1),pattern(iCell,:,i)));
            comp_patt_corr(i,iCell) = triu2vec(corrcoef(component(iCell,:,i),pattern(iCell,:,i)));
        end
    end
    
    Rp = ((patt_corr)-(comp_corr.*comp_patt_corr))./sqrt((1-comp_corr.^2).*(1-comp_patt_corr.^2));
    Rc = ((comp_corr)-(patt_corr.*comp_patt_corr))./sqrt((1-patt_corr.^2).*(1-comp_patt_corr.^2));
    Zp = (0.5.*log((1+Rp)./(1-Rp)))./sqrt(1./(nTestDirs-3));
    Zc = (0.5.*log((1+Rc)./(1-Rc)))./sqrt(1./(nTestDirs-3));
    
    [maxVal maxDir] = max(avg_resp_dir(:,:,1,1,1),[],2);
    oppDir = maxDir+nTestDirs/2;
    oppDir(find(oppDir>nTestDirs)) = oppDir(find(oppDir>nTestDirs))-nTestDirs;
    oppVal = indOnly(avg_resp_dir(:,:,1,1,1),oppDir);
    oppVal(find(oppVal<0)) = 0;
    DSI = (maxVal-oppVal)./(maxVal+oppVal);
    
    Zp_all = [Zp_all Zp];
    Zc_all = [Zc_all Zc];
        
    avg_resp_dir_shift = avg_resp_dir;
    avg_resp_dir_shift(:,:,2,:,:) = circshift(avg_resp_dir_shift(:,:,2,:,:),(maskDiffs/2)./int,2);
    avg_resp_dir_shift_all = cat(1, avg_resp_dir_shift_all, avg_resp_dir_shift);

    maxDir_all = [maxDir_all; maxDir];
    maxVal_all = [maxVal_all; maxVal];
    DSI_all = [DSI_all; DSI];

    h_resp_all = cat(1,h_resp_all,h_resp);
end

resp_ind_dir = find(sum(h_resp_all(:,:,1),2));
resp_ind_dir_use = intersect(resp_ind_dir,find(DSI_all>0.5));
    
Zc_use = intersect(resp_ind_dir, find(Zc_all(1,:)>1.28 & Zc_all(1,:)-Zp_all(1,:)>1.28));
Zp_use = intersect(resp_ind_dir, find(Zp_all(1,:)>1.28 & Zp_all(1,:)-Zc_all(1,:)>1.28));

figure;
start = 0;
for i = 1:5
ind_use = intersect(resp_ind_dir_use, find(maxDir_all==i));
subplot(5,2,1+start)
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_use,:,1,1,1),1) mean(avg_resp_dir_shift_all(ind_use,1,1,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_use,:,1,2,1),1) mean(avg_resp_dir_shift_all(ind_use,1,1,2,1),1)])
if start == 0
    title('Gratings')
end
subplot(5,2,2+start)
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_use,:,2,1,1),1) mean(avg_resp_dir_shift_all(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_use,:,2,2,1),1) mean(avg_resp_dir_shift_all(ind_use,1,2,2,1),1)])
if start == 0
    title('Plaids')
end
start = start+2;
end

print(fullfile(base,saveLoc, 'Analysis\2P\CrossOri\RandDirAdaptSummary', 'randDirAdapt_polarByPref.pdf'),'-dpdf','-bestfit')


ind = [3 11];
ind_orth = [];
ind_orth_Zc = [];
ind_orth_Zp = [];
for i = ind
    ind_orth = [ind_orth; intersect(resp_ind_dir_use, find(maxDir_all==i))];
    ind_orth_Zc = [ind_orth_Zc; intersect(Zc_use,intersect(resp_ind_dir_use, find(maxDir_all==i)))];
    ind_orth_Zp = [ind_orth_Zp; intersect(Zp_use,intersect(resp_ind_dir_use, find(maxDir_all==i)))];
end

figure;
[n n2] = subplotn(length(ind_orth_Zc));
for i = 1:length(ind_orth_Zc)
    subplot(n,n2,i)
    polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_orth_Zc(i),:,1,1,1),1) mean(avg_resp_dir_shift_all(ind_orth_Zc(i),1,1,1,1),1)])
    hold on
    polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_orth_Zc(i),:,1,2,1),1) mean(avg_resp_dir_shift_all(ind_orth_Zc(i),1,1,2,1),1)])
end
suptitle('Zc cells - Gratings')
figure;
for i = 1:length(ind_orth_Zc)
    subplot(n,n2,i)
    polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_orth_Zc(i),:,2,1,1),1) mean(avg_resp_dir_shift_all(ind_orth_Zc(i),1,2,1,1),1)])
    hold on
    polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_orth_Zc(i),:,2,2,1),1) mean(avg_resp_dir_shift_all(ind_orth_Zc(i),1,2,2,1),1)])
end
suptitle('Zc cells - Plaids')

figure;
[n n2] = subplotn(length(ind_orth_Zp));
for i = 1:length(ind_orth_Zp)
    subplot(n,n2,i)
    polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_orth_Zp(i),:,1,1,1),1) mean(avg_resp_dir_shift_all(ind_orth_Zp(i),1,1,1,1),1)])
    hold on
    polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_orth_Zp(i),:,1,2,1),1) mean(avg_resp_dir_shift_all(ind_orth_Zp(i),1,1,2,1),1)])
end
suptitle('Zp cells - Gratings')
print(fullfile(base, saveLoc,'Analysis\2P\CrossOri\RandDirAdaptSummary', 'randDirAdapt_exZpGratings.pdf'),'-dpdf','-bestfit')
figure;
for i = 1:length(ind_orth_Zp)
    subplot(n,n2,i)
    polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_orth_Zp(i),:,2,1,1),1) mean(avg_resp_dir_shift_all(ind_orth_Zp(i),1,2,1,1),1)])
    hold on
    polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_shift_all(ind_orth_Zp(i),:,2,2,1),1) mean(avg_resp_dir_shift_all(ind_orth_Zp(i),1,2,2,1),1)])
end
suptitle('Zp cells - Plaids')
print(fullfile(base, saveLoc,'Analysis\2P\CrossOri\RandDirAdaptSummary', 'randDirAdapt_exZpPlaids.pdf'),'-dpdf','-bestfit')

avg_resp_dir_align_all = avg_resp_dir_shift_all;
totCells = size(avg_resp_dir_align_all,1);
for iCell = 1:totCells
    avg_resp_dir_align_all(iCell,:,:,:,:) = circshift(avg_resp_dir_shift_all(iCell,:,:,:,:),4-maxDir_all(iCell),2);
end

tests = [12 1 2 6 7 8];
orth = [3 4 5 9 10 11];
tests_ind = [];
for i = 1:length(tests)
    tests_ind = [tests_ind; find(maxDir_all == tests(i))];
end
orth_ind = [];
for i = 1:length(orth)
    orth_ind = [orth_ind; find(maxDir_all == orth(i))];
end


ind_use = intersect(resp_ind_dir_use, [find(maxDir_all==11); find(maxDir_all==3)]);
figure;
subplot(2,2,1)
scatter(Zc_all(1,ind_use), Zp_all(1,ind_use))
ylim([-4 8])
xlim([-4 8])
axis square
hold on
plotZcZpBorders
xlabel('Zc')
ylabel('Zp')
subplot(2,2,2)
scatter(Zc_all(2,ind_use), Zp_all(2,ind_use))
ylim([-4 8])
xlim([-4 8])
axis square
hold on
plotZcZpBorders
xlabel('Zc')
ylabel('Zp')
subplot(2,2,3)
scatter(Zc_all(1,ind_use), Zc_all(2,ind_use))
ylim([-4 8])
xlim([-4 8])
axis square
refline(1,0)
[h p] = ttest(Zc_all(1,ind_use), Zc_all(2,ind_use));
xlabel('Zc- control')
ylabel('Zc- adapt')
title(['p = ' num2str(chop(p,2))])
subplot(2,2,4)
scatter(Zp_all(1,ind_use), Zp_all(2,ind_use))
ylim([-4 8])
xlim([-4 8])
axis square
refline(1,0)
[h p] = ttest(Zp_all(1,ind_use), Zp_all(2,ind_use));
xlabel('Zp- control')
ylabel('Zp- adapt')
title(['p = ' num2str(chop(p,2))])
suptitle(['Pref: 60 & 300- n = ' num2str(length(ind_use))])
print(fullfile(base, saveLoc, 'Analysis\2P\CrossOri\RandDirAdaptSummary', 'randDirAdapt_ZcZpEffects.pdf'),'-dpdf','-bestfit')

ind_use = intersect(resp_ind_dir_use, orth_ind);
figure;
subplot(2,2,1)
scatter(Zc_all(1,ind_use), Zp_all(1,ind_use))
ylim([-4 8])
xlim([-4 8])
axis square
hold on
plotZcZpBorders
xlabel('Zc')
ylabel('Zp')
subplot(2,2,2)
scatter(Zc_all(2,ind_use), Zp_all(2,ind_use))
ylim([-4 8])
xlim([-4 8])
axis square
hold on
plotZcZpBorders
xlabel('Zc')
ylabel('Zp')
subplot(2,2,3)
scatter(Zc_all(1,ind_use), Zc_all(2,ind_use))
ylim([-4 8])
xlim([-4 8])
axis square
refline(1,0)
[h p] = ttest(Zc_all(1,ind_use), Zc_all(2,ind_use));
xlabel('Zc- control')
ylabel('Zc- adapt')
title(['p = ' num2str(chop(p,2))])
subplot(2,2,4)
scatter(Zp_all(1,ind_use), Zp_all(2,ind_use))
ylim([-4 8])
xlim([-4 8])
axis square
refline(1,0)
[h p] = ttest(Zp_all(1,ind_use), Zp_all(2,ind_use));
xlabel('Zp- control')
ylabel('Zp- adapt')
title(['p = ' num2str(chop(p,2))])
suptitle(['Pref: 60, 90 120, 240, 270, 300- n = ' num2str(length(ind_use))])
print(fullfile(base, saveLoc, 'Analysis\2P\CrossOri\RandDirAdaptSummary', 'randDirAdapt_ZcZpEffects_allOrth.pdf'),'-dpdf','-bestfit')


figure;
subplot(3,2,1)
ind_use = intersect(resp_ind_dir_use, tests_ind);
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,1,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,1,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,1,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,1,2,1),1)])
title(num2str(length(ind_use)))
ind_use = intersect(Zc_use,intersect(resp_ind_dir_use, tests_ind));
subplot(3,2,3)
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,2,1),1)])
title(num2str(length(ind_use)))
ind_use = intersect(Zp_use,intersect(resp_ind_dir_use, tests_ind));
subplot(3,2,5)
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,2,1),1)])
title(num2str(length(ind_use)))
subplot(3,2,2)
ind_use = intersect(resp_ind_dir_use, orth_ind);
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,1,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,1,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,1,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,1,2,1),1)])
title(num2str(length(ind_use)))
ind_use = intersect(Zc_use,intersect(resp_ind_dir_use, orth_ind));
subplot(3,2,4)
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,2,1),1)])
title(num2str(length(ind_use)))
ind_use = intersect(Zp_use,intersect(resp_ind_dir_use, orth_ind));
subplot(3,2,6)
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,2,1),1)])
title(num2str(length(ind_use)))
suptitle({'Left: 0, 30, 150, 180, 210, 330; Right: 60, 90, 120, 240, 270, 300;' 'Top: Grating; Middle: Zc; Bottom: Zp'})
print(fullfile(base, saveLoc, 'Analysis\2P\CrossOri\RandDirAdaptSummary', 'randDirAdapt_Polar_adaptVortho.pdf'),'-dpdf','-bestfit')

figure;
subplot(3,2,1)
ind_use = intersect(resp_ind_dir_use, [find(maxDir_all==1); find(maxDir_all==7)]);
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,1,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,1,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,1,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,1,2,1),1)])
title(num2str(length(ind_use)))
ind_use = intersect(Zc_use,intersect(resp_ind_dir_use, [find(maxDir_all==1); find(maxDir_all==7)]));
subplot(3,2,3)
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,2,1),1)])
title(num2str(length(ind_use)))
ind_use = intersect(Zp_use,intersect(resp_ind_dir_use, [find(maxDir_all==1); find(maxDir_all==7)]));
subplot(3,2,5)
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,2,1),1)])
title(num2str(length(ind_use)))
subplot(3,2,2)
ind_use = intersect(resp_ind_dir_use, [find(maxDir_all==11); find(maxDir_all==3)]);
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,1,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,1,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,1,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,1,2,1),1)])
title(num2str(length(ind_use)))
ind_use = intersect(Zc_use,intersect(resp_ind_dir_use, [find(maxDir_all==11); find(maxDir_all==3)]));
subplot(3,2,4)
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,2,1),1)])
title(num2str(length(ind_use)))
ind_use = intersect(Zp_use,intersect(resp_ind_dir_use, [find(maxDir_all==11); find(maxDir_all==3)]));
subplot(3,2,6)
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,1,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,1,1),1)])
hold on
polar(deg2rad([testDirs testDirs(1)]), [mean(avg_resp_dir_align_all(ind_use,:,2,2,1),1) mean(avg_resp_dir_align_all(ind_use,1,2,2,1),1)])
title(num2str(length(ind_use)))
suptitle({'Left: 0 & 180; Right: 60 & 300;' 'Top: Grating; Middle: Zc; Bottom: Zp'})
print(fullfile(base, saveLoc, 'Analysis\2P\CrossOri\RandDirAdaptSummary', 'randDirAdapt_Polar_0v60.pdf'),'-dpdf','-bestfit')