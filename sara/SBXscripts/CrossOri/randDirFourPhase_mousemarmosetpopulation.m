clc;clear all;close all
%% pull mouse data
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC'); 
area = 'all_areas';
area_list = strvcat('V1');
narea = length(area_list);
nCells = [];


for iA = 1
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];

    pattern_ind = (Zp_all-Zc_all);    
    pattind = mean(pattern_ind,1);
    pattpeak = max(pattern_ind,[],1);
    
    figure(1)
        subplot(4,4,1)
            scatter(pattpeak(resp_ind),amp_all(resp_ind),10,'filled')
            hold on
            xlabel('Peak pattern index (Zp-Zc)')
            ylabel('Spatial variance (amp)')
            xlim([-4 8])
            set(gca,'TickDir','out'); axis square
            subtitle(leg_str(:,iA))
end

%% load marm data
clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';


load(fullfile(base, (['Data\fromNicholas\lindsey_frg1.mat'])))

x=[-150:30:180];
x_rad = deg2rad(x);

stimDirs = [0 30 60 90 120 150 180 210 240 270 300 330];
nStimDir = 12;
nMaskPhas = 2;
nCells = 3;



%Fit PCI with sinusoid

PCI = (Zpd-Zcd);
nCells = size(PCI,2);

%fit sinusoid
phase = [0 90 180 270];
phase_range = 0:1:359;

for iCell = 1:nCells
    [b_hat_all(iCell,1), amp_hat_all(iCell,1), per_hat_all(iCell,1),pha_hat_all(iCell,1),sse_all(iCell,1),R_square_all(iCell,1)] = sinefit_PCI(deg2rad(phase),PCI(:,iCell));
    yfit_all(iCell,:,1) = b_hat_all(iCell,1)+amp_hat_all(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(iCell,1) + 2.*pi/pha_hat_all(iCell,1)));
end


%% Plot population-- PCI by variance (amplitude)

pattern_ind = (Zpd-Zcd);    
pattind = mean(pattern_ind,1);
pattpeak = max(pattern_ind,[],1);

figure(1)
    subplot(4,4,1)
    scatter(pattpeak,amp_hat_all,10,'filled')
    hold on



%% print figs
    figure(1); print(fullfile(outDir, [svName '_mousemarmoset_population.pdf']),'-dpdf', '-fillpage') 
