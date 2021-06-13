clear all;
close all;
clc;
doRedChannel = 0;
ds = 'CrossOriRandDir_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'PhaseRev');
phaseRev = struct;
area_list = {'V1'};
driver_list = {'SLC','SOM'};
narea =length(area_list);
ndriver =length(driver_list);
for id = 1:ndriver
    driver = driver_list{id};
    area = area_list{1};
    fprintf([driver '\n'])
    totCells = 0;
    stim_OSI_all = [];
    k1_all = [];
    resp_ind_all = [];
    f1_all = [];
    f2_all = [];
    f2overf1_all = [];
    mouse_list = [];

    for iexp = 1:nexp
        if sum(strcmp(expt(iexp).img_loc,area)) & ~isempty(expt(iexp).prFolder) & strcmp(expt(iexp).driver,driver)        
            mouse = expt(iexp).mouse;
            mouse_list = strvcat(mouse_list, mouse);
            date = expt(iexp).date;
            ImgFolder = expt(iexp).coFolder;
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);
            fprintf([mouse ' ' date '\n'])
            
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']))
            nCells = size(h_resp,1);
            if length(expt(iexp).img_loc)>1
                i = find(strcmp(expt(iexp).img_loc,area));
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_splitImage.mat']))
                ind = find(maskCat==i);
                h_resp = h_resp(ind,:,:);
                nCells = length(ind);
                stim_OSI = stim_OSI(1,ind);
                k1_dir = k1_dir(1,ind);
            end
            stim_OSI_all = [stim_OSI_all stim_OSI];
            k1_all = [k1_all k1_dir];
            resp_ind = find(sum(sum(h_resp,2),3));
            resp_ind_all = [resp_ind_all resp_ind'+totCells];
            
            ImgFolder = expt(iexp).prFolder;
            nrun = length(ImgFolder);
            run_str = catRunName(cell2mat(ImgFolder), nrun);
            load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_f1f2.mat']))
            if length(expt(iexp).img_loc)>1
                f1 = f1(ind);
                f2 = f2(ind);
                f2overf1 = f2overf1(ind);
            end
            f1_all = [f1_all f1];
            f2_all = [f2_all f2];
            f2overf1_all = [f2overf1_all f2overf1];
            totCells = totCells+nCells;
        end
    end
    phaseRev(id).name = area;
    phaseRev(id).nExpt = size(mouse_list,1);
    phaseRev(id).nMice = size(unique(mouse_list,'rows'),1);
    phaseRev(id).totCells = totCells;
    phaseRev(id).respCells = length(resp_ind_all);
    phaseRev(id).threshCells = length(find(f1_all>0.02));
    
    save(fullfile(summaryDir,['phaseRev_Summary_' driver '.mat']),'resp_ind_all', 'f1_all','f2_all','f2overf1_all','mouse_list')
    figure(1)
    subplot(2,3,1)
    cdfplot(f1_all(find(f1_all>0.02)))
    hold on
    subplot(2,3,2)
    cdfplot(f2_all(find(f1_all>0.02)))
    hold on
    subplot(2,3,3)
    cdfplot(f2overf1_all(find(f1_all>0.02)))
    hold on
    subplot(2,2,4)
    errorbar(id, mean(f2overf1_all(find(f1_all>0.02))),std(f2overf1_all(find(f1_all>0.02)))./sqrt(phaseRev(id).threshCells),'o');
    hold on
    subplot(2,2,3)
    errorbar(mean(f1_all(find(f1_all>0.02))),mean(f2_all(find(f1_all>0.02))),std(f1_all(find(f1_all>0.02)))./sqrt(phaseRev(id).threshCells),std(f1_all(find(f1_all>0.02)))./sqrt(phaseRev(id).threshCells),std(f2_all(find(f1_all>0.02)))./sqrt(phaseRev(id).threshCells),std(f2_all(find(f1_all>0.02)))./sqrt(phaseRev(id).threshCells),'o');
    hold on
    leg{id} = [driver '- n = ' num2str(phaseRev(id).threshCells)];
    
    figure(2)
    subplot(2,ndriver,id)
    ind_h = intersect(find(f1_all>0.02),find(stim_OSI_all>0.5));
    ind_l = intersect(find(f1_all>0.02),find(stim_OSI_all<=0.5));
    cdfplot(f2overf1_all(ind_h))
    hold on
    cdfplot(f2overf1_all(ind_l))
    legend({['OSI>0.5 - n = ' num2str(length(ind_h))],['OSI<0.5 - n = ' num2str(length(ind_l))]}, 'location', 'southeast')
    xlabel('F2/F1')
    title(driver)
    xlim([0 2.5])
    subplot(2,ndriver,id+ndriver)
    ind_h = intersect(find(f1_all>0.02),find(k1_all>12));
    ind_l = intersect(find(f1_all>0.02),find(k1_all<=12));
    cdfplot(f2overf1_all(ind_h))
    hold on
    cdfplot(f2overf1_all(ind_l))
    legend({['K>12 - n = ' num2str(length(ind_h))],['K<12 - n = ' num2str(length(ind_l))]}, 'location', 'southeast')
    xlabel('F2/F1')
    xlim([0 2.5])
    title(driver)
    
    figure(3)
    subplot(2,ndriver,id)
    ind_h = intersect(find(f1_all>0.02),find(f2overf1_all>0.5));
    ind_l = intersect(find(f1_all>0.02),find(f2overf1_all<=0.5));
    cdfplot(stim_OSI_all(ind_h))
    hold on
    cdfplot(stim_OSI_all(ind_l))
    legend({['F2/F1>0.5 - n = ' num2str(length(ind_h))],['F2/F1<0.5 - n = ' num2str(length(ind_l))]}, 'location', 'northwest')
    xlabel('OSI')
    title(driver)
    xlim([0 1])
    subplot(2,ndriver,id+ndriver)
    ind_h = intersect(find(f1_all>0.02),find(f2overf1_all>0.5));
    ind_l = intersect(find(f1_all>0.02),find(f2overf1_all<=0.5));
    cdfplot(k1_all(ind_h))
    hold on
    cdfplot(k1_all(ind_l))
    legend({['F2/F1>0.5 - n = ' num2str(length(ind_h))],['F2/F1<0.5 - n = ' num2str(length(ind_l))]}, 'location', 'northwest')
    xlabel('Tuning Sharpness (k)')
    xlim([0 30])
    title(driver)
end

phaseRevTable = struct2table(phaseRev);
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
stim.frameRate = frame_rate;
stim.nOnSec = double(input.nScansOn)./stim.frameRate;
stim.nOffSec = double(input.nScansOff)./stim.frameRate;
stim.phaseRevHz = stim.frameRate./double(input.nScansPhaseCyc.*2);
stim.gratingDiameter = input.gratingDiameterDeg;
stim.gratingSF = input.gratingSpatialFreqCPD;
stim.gratingOriStep = input.gratingDirectionStepDeg;
stim.gratingOriStepN = input.gratingDirectionStepN;
stim.gratingPhaseStep = input.gratingStartingPhaseStepDeg;
stim.gratingPhaseStepN = input.gratingStartingPhaseStepN;

save(fullfile(summaryDir,['phaseRevStats_allDriver.mat']),'stim','phaseRev');

figure(1)
subplot(2,3,1)
xlabel('F1')
ylabel('Fraction of cells')
legend(leg,'location','southeast')
xlim([0 0.5])
vline(0.02)
title('')
subplot(2,3,2)
xlabel('F2')
ylabel('Fraction of cells')
legend(leg,'location','southeast')
xlim([0 0.5])
title('')
subplot(2,3,3)
xlabel('F2/F1')
ylabel('Fraction of cells')
legend(leg,'location','southeast')
xlim([0 2])
title('')
subplot(2,2,4)
xlim([0 ndriver+1])
ylim([0 1])
set(gca,'XTick',1:ndriver,'XTickLabel',driver_list)
ylabel('F2/F1')
subplot(2,2,3)
xlabel('F1')
ylabel('F2')
xlim([0 0.15])
ylim([0 0.15])
print(fullfile(summaryDir,['phaseRev_Driver_Summary.pdf']),'-dpdf','-bestfit')

figure(2)
print(fullfile(summaryDir,['phaseRev_Driver_TuningVsF2F1.pdf']),'-dpdf','-bestfit')
figure(3)
print(fullfile(summaryDir,['phaseRev_Driver_TuningVsF2F1_invert.pdf']),'-dpdf','-bestfit')
