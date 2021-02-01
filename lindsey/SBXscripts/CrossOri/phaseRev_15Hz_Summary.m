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
area_list = ['V1'; 'LM'; 'AL'; 'RL'; 'PM'];
driver = 'SLC';
narea =length(area_list);
figure;
for iarea = 1:narea
    area = area_list(iarea,:);
    fprintf([area '\n'])
    totCells = 0;
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
            nCells = size(h_resp,1);
            if length(expt(iexp).img_loc)>1
                i = find(strcmp(expt(iexp).img_loc,area));
                load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_splitImage.mat']))
                ind = find(maskCat==i);
                h_resp = h_resp(ind,:,:);
                nCells = length(ind);
            end
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
    phaseRev(iarea).name = area;
    phaseRev(iarea).nExpt = size(mouse_list,1);
    phaseRev(iarea).nMice = size(unique(mouse_list,'rows'),1);
    phaseRev(iarea).totCells = totCells;
    phaseRev(iarea).respCells = length(resp_ind_all);
    phaseRev(iarea).threshCells = length(intersect(resp_ind_all,find(f1_all>0.02)));
    
    save(fullfile(summaryDir,['phaseRev_Summary_' area '.mat']),'resp_ind_all', 'f1_all','f2_all','f2overf1_all','mouse_list')
    subplot(3,2,1)
    cdfplot(f2overf1_all(intersect(resp_ind_all,find(f1_all>0.02))))
    hold on
    subplot(3,2,3)
    errorbar(iarea, mean(f2overf1_all(intersect(resp_ind_all,find(f1_all>0.02)))),std(f2overf1_all(intersect(resp_ind_all,find(f1_all>0.02))))./sqrt(phaseRev(iarea).threshCells),'o');
    hold on
    subplot(3,2,5)
    errorbar(mean(f1_all(intersect(resp_ind_all,find(f1_all>0.02)))),mean(f2_all(intersect(resp_ind_all,find(f1_all>0.02)))),std(f1_all(intersect(resp_ind_all,find(f1_all>0.02))))./sqrt(phaseRev(iarea).threshCells),std(f1_all(intersect(resp_ind_all,find(f1_all>0.02))))./sqrt(phaseRev(iarea).threshCells),std(f2_all(intersect(resp_ind_all,find(f1_all>0.02))))./sqrt(phaseRev(iarea).threshCells),std(f2_all(intersect(resp_ind_all,find(f1_all>0.02))))./sqrt(phaseRev(iarea).threshCells),'o');
    hold on
    leg{iarea} = [area '- n = ' num2str(phaseRev(iarea).threshCells)];
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

save(fullfile(summaryDir,['phaseRevStats_allArea.mat']),'stim','phaseRev');

subplot(3,2,1)
xlabel('F2/F1')
ylabel('Fraction of cells')
legend(leg,'location','southeast')
xlim([0 2])
title('15Hz')
subplot(3,2,3)
xlim([0 narea+1])
ylim([0 1])
set(gca,'XTick',1:narea,'XTickLabel',area_list)
ylabel('F2/F1')
subplot(3,2,5)
xlabel('F1')
ylabel('F2')
xlim([0 0.15])
ylim([0 0.05])

load(fullfile(summaryDir,['phaseRevStats_V1only_30Hz.mat']))
load(fullfile(summaryDir,['phaseRev30Hz_Summary.mat']))
subplot(3,2,2)
cdfplot(f2overf1_all(intersect(resp_ind_all,find(f1_all>0.02))))
hold on
xlabel('F2/F1')
ylabel('Fraction of cells')
legend(['V1- n = ' num2str(phaseRev.threshCells)],'location','southeast')
xlim([0 2])
title('30 Hz')
subplot(3,2,4)
errorbar(1, mean(f2overf1_all(intersect(resp_ind_all,find(f1_all>0.02)))),std(f2overf1_all(intersect(resp_ind_all,find(f1_all>0.02))))./sqrt(phaseRev.threshCells),'o');
hold on
xlim([0 2])
ylim([0 1])
set(gca,'XTick',1,'XTickLabel','V1')
ylabel('F2/F1')
subplot(3,2,6)
errorbar(mean(f1_all(intersect(resp_ind_all,find(f1_all>0.02)))),mean(f2_all(intersect(resp_ind_all,find(f1_all>0.02)))),std(f1_all(intersect(resp_ind_all,find(f1_all>0.02))))./sqrt(phaseRev.threshCells),std(f1_all(intersect(resp_ind_all,find(f1_all>0.02))))./sqrt(phaseRev.threshCells),std(f2_all(intersect(resp_ind_all,find(f1_all>0.02))))./sqrt(phaseRev.threshCells),std(f2_all(intersect(resp_ind_all,find(f1_all>0.02))))./sqrt(phaseRev.threshCells),'o');
hold on
xlabel('F1')
ylabel('F2')
xlim([0 0.15])
ylim([0 0.05])

print(fullfile(summaryDir,['phaseRev_Summary.pdf']),'-dpdf','-bestfit')
