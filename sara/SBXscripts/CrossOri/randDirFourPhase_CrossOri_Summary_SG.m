close all; clear all; clc;
doRedChannel = 0;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');

doPlot = 1;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
svName = 'randPhase';
eval(ds)
driver = 'SLC'; %Scn
img_area = {'PM';'L4'}; %LM
inj_area = 'PM';

max_dist = 10;

rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);

sig_stim = [];
sig_dir = [];
prefDir_resamp_all = [];

Zp_all = [];
Zc_all = [];
PCI_all = [];
ZpZcPWdist_all = [];

DSI_all = [];

yfit_all_all = [];
b_all = [];
amp_all = [];
pha_all = [];
per_all = [];
sse_all_all = [];
Rsq_all  = [];

b_all_sh = [];
amp_all_sh = [];
Rsq_all_sh  = [];
sse_all_all_sh = [];

dir_b_all = [];
k1_all = [];
R1_all = [];
R2_all = [];
u1_all = [];
u2_all = [];
dir_Rsq_all = [];
dir_sse_all = [];
dir_yfits = [];

plaid_corr_all = [];

std_avg_all = [];
mean_avg_all = [];
grat_pref_all = [];
plaid_pref_all = [];

g_dsi_all = [];
g_osi_all = [];
ang_dir_all = [];
ang_ori_all = [];

red_cells_all = [];

mouse_list = [];
totCells = [];

% V1 L2/3 -- 9 10 13 26 40 51 53 67 75
% V1 L4 -- 63 64 107 109 113 114 115 116
% AL L2/3 -- 17 27 58 70 74 90 92
% LM L2/3 -- 19 30 60 73 78 87 88 91
% PM L2/3 -- 97 102 103 104 122

start=1;
for iexp = [97 102 103 104 122] 
    mouse = expt(iexp).mouse;
    mouse_list = strvcat(mouse_list, mouse);
    date = expt(iexp).date;
    if isfield(expt,'copFolder') 
        ImgFolder = expt(iexp).copFolder;
    else
        ImgFolder = expt(iexp).coFolder;
    end
    nrun = length(ImgFolder);
    run_str = 'runs-002';
        
    fprintf([mouse ' ' date '\n'])
        
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_plaidCorr.mat']))
    if exist(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits.mat']))
        load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits.mat']))
    else
        fprintf(['PCI fit file does not exist for    ' date '_' mouse '  max_dist ' num2str(max_dist) '\n'])
        stop
    end
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str],['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_ZpZc_pairwiseDist.mat']))
    
    if doRedChannel==1
        load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']),'red_cells')
    else
        red_cells=0;
    end
    
    nCells = size(resp_cell{end,end,end},1);
    
    if start>1
        red_cells_all = [red_cells_all; red_cells+totCells];
    else
        red_cells_all = red_cells;
    end
    
    plaid_corr_all = [plaid_corr_all; plaid_corr'];
    
    b_all = [b_all; b_hat_all];
    amp_all = [amp_all; amp_hat_all];
    Rsq_all = [Rsq_all; R_square_all];
    sse_all_all = [sse_all_all; sse_all];
    
    b_all_sh = [b_all_sh; b_hat_all];
    amp_all_sh = [amp_all_sh; amp_hat_all];
    Rsq_all_sh = [Rsq_all_sh; R_square_all];
    sse_all_all_sh = [sse_all_all_sh; sse_all];
    
    yfit_all_all = cat(1,yfit_all_all, yfit_all);
    
    if start > 1
        sig_stim = [sig_stim; resp_ind_dir + totCells]; 
        sig_dir = [sig_dir, p_dir + totCells];   
    else
        sig_stim = resp_ind_dir; 
        sig_dir = p_dir;   
    end
    
    prefDir_resamp_all = [prefDir_resamp_all; prefDir_resamp'];
    
    if start > 1
        totCells = totCells + nCells;    
    else
        totCells = nCells;
    end
    
    Zp_all = cat(2,Zp_all, Zp);
    Zc_all = cat(2,Zc_all, Zc);
    PCI_all = cat(2,PCI_all, PCI);
    ZpZcPWdist_all = cat(2,ZpZcPWdist_all, ZpZcPWdist);
    
    DSI_all = [DSI_all; DSI'];

    nStimDir=size(resp_cell,1);
    nMaskPhas=size(resp_cell,2);
    std_bystim=zeros(nCells,nStimDir,nMaskPhas+1);
    mean_bystim=zeros(nCells,nStimDir,nMaskPhas+1);
    for iDir = 1:nStimDir
        mean_bystim(:,iDir,1) = mean(resp_cell{iDir,1,1},2);
        std_bystim(:,iDir,1) = std(resp_cell{iDir,1,1},0,2);
        for ip = 1:nMaskPhas
            mean_bystim(:,iDir,ip+1) = mean(resp_cell{iDir,ip,2},2);
            std_bystim(:,iDir,ip+1) = std(resp_cell{iDir,ip,2},0,2);
        end
    end
    mean_avg = mean(mean(mean(mean_bystim,3),4),2);
    std_avg = mean(mean(mean(std_bystim,3),4),2);
    std_avg_all = [std_avg_all; std_avg];
    mean_avg_all = [mean_avg_all; mean_avg];
    grat_resp = mean_bystim(:,1);
    plaid_resp = mean(mean(mean_bystim(:,2:5),2),3);
    grat_pref = find(grat_resp>plaid_resp);
    plaid_pref = find(plaid_resp>grat_resp);
    grat_pref_all = [grat_pref_all; grat_pref+ totCells];
    plaid_pref_all = [plaid_pref_all; plaid_pref+ totCells];


    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], ['centroid_dist_max_' num2str(max_dist)], [date '_' mouse '_' run_str '_DirectionTuningFit.mat']))

    dir_b_all = [dir_b_all; b_hat_all];
    k1_all = [k1_all; k1_hat_all];
    R1_all = [R1_all; R1_hat_all];
    R2_all = [R2_all; R2_hat_all];
    u1_all = [u1_all; u1_hat_all];
    u2_all = [u2_all; u2_hat_all];
    dir_Rsq_all = [dir_Rsq_all; R_square_all];
    dir_sse_all = [dir_sse_all; sse_all];
    % dir_yfits = cat(1,dir_yfits, dir_yfit_all');
    
    cellnum(start,1) = nCells;
    cellnum(start,2) = totCells;
    cellnum(start,3) = length(amp_hat_all);

    g_dsi_all = [g_dsi_all; g_dsi'];
    g_osi_all = [g_osi_all; g_osi'];
    ang_dir_all = [ang_dir_all; ang_dir'];
    ang_ori_all = [ang_ori_all; ang_ori'];
    
    start=start+1;
end
    save(fullfile(summaryDir,[svName '_Summary_' inj_area  '_' driver '.mat']), 'ang_dir_all', 'ang_ori_all', 'g_dsi_all', 'g_osi_all', 'grat_pref_all', 'plaid_pref_all', 'mean_avg_all', 'std_avg_all', 'prefDir_resamp_all', 'sig_dir','sig_stim','plaid_corr_all', 'red_cells_all', 'ZpZcPWdist_all', 'b_all_sh', 'amp_all_sh', 'Rsq_all_sh', 'sse_all_all_sh','dir_yfits','R2_all','u2_all','dir_b_all', 'k1_all', 'R1_all', 'u1_all', 'dir_Rsq_all', 'dir_sse_all', 'nCells', 'Zp_all', 'Zc_all', 'PCI_all', 'b_all','amp_all','Rsq_all','yfit_all_all', 'sse_all_all', 'DSI_all','mouse_list')
    %save(fullfile(summaryDir,[svName '_Summary_' inj_area '_' num2str(iexp) '.mat']), 'ang_dir_all', 'ang_ori_all', 'g_dsi_all', 'g_osi_all', 'mean_avg_all', 'std_avg_all', 'prefDir_resamp_all','sig_dir','sig_stim','plaid_corr_all', 'red_cells_all', 'dir_yfits','ZpZcPWdist_all','R2_all','u2_all','dir_b_all', 'k1_all', 'R1_all', 'u1_all', 'dir_Rsq_all', 'dir_sse_all', 'nCells', 'Zp_all', 'Zc_all', 'PCI_all', 'b_all','amp_all','Rsq_all','yfit_all_all', 'sse_all_all', 'DSI_all','mouse_list')

%%
close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary','summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC','SLC','SLC'); 
area = 'all_areas';
area_list = strvcat('V1', 'LM', 'AL');
narea = length(area_list);
nCells = [];


figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];
    

    subplot(3,2,1)
        cdfplot(Rsq_all(resp_ind,:))
        hold on
        xlabel('fit Rsq')
        legend(leg_str, 'location', 'southeast'); set(gca,'TickDir','out'); box off; grid off
        
    subplot(3,2,2)
        cdfplot(sse_all_all(resp_ind,:))
        hold on
        xlabel('fit SSE')
        legend(leg_str, 'location', 'southeast'); set(gca,'TickDir','out'); box off; grid off
    
    subplot(3,2,3)
        cdfplot(b_all(resp_ind,:))
        hold on
        xlabel('fit baseline')
        legend(leg_str, 'location', 'southeast'); set(gca,'TickDir','out'); box off; grid off
    
    subplot(3,2,4)
        cdfplot(amp_all(resp_ind,:))
        hold on
        xlabel('fit amplitude')
        legend(leg_str, 'location', 'southeast'); set(gca,'TickDir','out'); box off; grid off
        
    subplot(3,2,5)
        cdfplot(reshape(Zp_all(:,resp_ind),[],1))
        hold on
        xlabel('Zp all')
        legend(leg_str, 'location', 'southeast'); set(gca,'TickDir','out'); box off; grid off
        
    subplot(3,2,6)
        cdfplot(reshape(PCI_all(:,resp_ind),[],1))
        hold on
        xlabel('PCI all')
        legend(leg_str, 'location', 'southeast'); set(gca,'TickDir','out'); box off; grid off
end

figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];
    for ip = 1:4
        subplot(narea,4,(narea*4-4)+ip)
            scatter(Zc_all(ip,resp_ind),Zp_all(ip,resp_ind))
            hold on
            plotZcZpBorders
            subtitle([area_list(iA,:)]); xlabel('Zc'); ylabel('Zp')
            set(gca,'TickDir','out'); box off; axis square; grid off
    end
end

    
    
%%
clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandDirFourPhase_ExptList_SG';
rc = behavConstsAV;
eval(ds)
nexp = length(expt);
frame_rate = 15;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randDir_randPhase';

figure;
for iexp = [1 2 4]
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc;
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    fprintf([mouse ' ' date '\n'])
    load(fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PatternCorrelationIndexFits.mat']))
    
    subplot(2,3,1)
        cdfplot(sse_all)
        hold on
        xlabel('summed square of residuals (SSE)')
        set(gca,'TickDir','out'); box off; axis square; grid off
        legend('exp1','exp2','exp4')
        
    subplot(2,3,2)
        cdfplot(b_hat_all)
        hold on
        xlabel('PCI (Zp-Zc) baseline')
        set(gca,'TickDir','out'); box off; axis square; grid off
        
    subplot(2,3,3)
        cdfplot(amp_hat_all)
        hold on
        xlabel('PCI (Zp-Zc) amplitude)')
        set(gca,'TickDir','out'); box off; axis square; grid off
        
    subplot(2,3,4)
        cdfplot(sse_all(ind))
        hold on
        xlabel('summed square of residuals (SSE)')
        subtitle('ind')
        set(gca,'TickDir','out'); box off; axis square; grid off
        
    subplot(2,3,5)
        cdfplot(b_hat_all(ind))
        hold on
        xlabel('PCI (Zp-Zc) baseline')
        subtitle('ind')
        set(gca,'TickDir','out'); box off; axis square; grid off
        
    subplot(2,3,6)
        cdfplot(amp_hat_all(ind))
        hold on
        xlabel('PCI (Zp-Zc) amplitude)')
        subtitle('ind')
        set(gca,'TickDir','out'); box off; axis square; grid off
end

print(fullfile(outDir, [svName '_PCIfitsbyExperiment.pdf']),'-dpdf', '-fillpage') 


