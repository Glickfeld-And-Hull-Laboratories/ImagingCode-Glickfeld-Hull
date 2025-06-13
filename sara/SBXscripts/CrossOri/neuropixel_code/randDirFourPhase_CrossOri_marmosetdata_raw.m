%% Analyze raw marmoset data sent from Nicholas Priebe, 4/18/2025
clc; clear all; close all;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
svName = 'randDirFourPhase_CrossOri';
bnsz = 10; %spike binsize = 10 ms

doPlot=1;
doGratPlot=0;

expts = strvcat('g01b', 'g06b', 'g12b', 'g17b'); %Four acute penetrations in V1, 1 marmoset
nexp = length(expts);

iexp = 3;

    fprintf([expts(iexp,:) '\n'])

    % Load data from experiment. Each experiment file contains one
    % variable, resp, that is nCells by nDir by nPhase by Grating/Plaid by
    % nTrials by Time (in 10 ms bins). For example, 310x12x4x2x13x120.
    % Grating type only has the first phase filled in.
    load(fullfile(base, 'Data\fromNicholas\CrossOri_randDirFourPhase_V1_marmoset', ([expts(iexp,:) '.mat'])))
    
    nCells = size(resp,1);
    nDirs = size(resp,2);
    nPhas = size(resp,3);
    nTrials = size(resp,5);
    nStim = (size(resp,3)+1)*size(resp,2);

    resp_cell = resp(:,:,:,:,:,21:120);
    base_cell = resp(:,:,:,:,:,1:20);
    avg_resp_dir(:,:,:,:,1) = mean(sum(resp_cell,6)/100,5); % Average across time and trial
    avg_resp_dir(:,:,:,:,2) = std(sum(resp_cell,6)/100,0,5)/nTrials; % SEM for each stim condition
    resp_dir_tc = mean(resp,5); % Average across trial, to look at time 
    resp_dir_tr = mean(resp,6); % Average across time, to look at trial 


    % Find cells that are significantly responsive to 1 stimulus
    resp_cell_trials(:,:,:,:,:,1) = sum(resp_cell,6)/100; % convert to one spike rate per trial (in Hz)
    base_cell_trials(:,:,:,:,:,1) = (sum(resp_cell,6)*5)/100; % multiply by 5 to convert to Hz
    
    for id = 1:nDirs
        [h_resp(:,id,1,1), p_resp(:,id,1,1)] = ttest2(squeeze(resp_cell_trials(:,id,1,1,:)),squeeze(base_cell_trials(:,id,1,1,:)),'dim',2,'tail','right','alpha', 0.05./nStim);
        for ip = 1:nPhas
            [h_resp(:,id,ip,2), p_resp(:,id,ip,2)] = ttest2(squeeze(resp_cell_trials(:,id,ip,2,:)),squeeze(base_cell_trials(:,id,ip,2,:)),'dim',2,'tail','right','alpha', 0.05./nStim);
        end
    end
    resp_ind_dir = find(sum(h_resp(:,:,1,1),2)); %sig responsive to gratings
    

    % figure; %trying to figure out if i have stim on + stim off or just stim on
    %     for ic = 1:36
    %         subplot(6,6,ic)
    %         plot(0:119,squeeze(resp_dir_tc(ic,1,1,1,:)))
    %     end


    % Run ANOVA across stim direction
    grat_resp = squeeze(resp_dir_tr(:,:,1,1,:)); % Make matrix nCells x nDir x nTrials
    plaid_resp = squeeze(resp_dir_tr(:,:,:,2,:)); % Make matrix nCells x nDir x nPhas x nTrials
    plaid_resp = permute(plaid_resp,[1,3,2,4]);

    group_mat = zeros(nDirs,nTrials);
    for id = 1:nDirs
        group_mat(id,:) = id;
    end

    all_grat_resp = reshape(grat_resp,[nCells,nDirs*nTrials]);
    all_plaid_resp = reshape(plaid_resp,[nCells,nPhas,nDirs*nTrials]);
    group = reshape(group_mat,[1,nDirs*nTrials]);

    p_anova_dir = zeros(1,nCells);
    for iCell = 1:nCells
        p_anova_dir(iCell) = anova1(all_grat_resp(iCell,:), group, 'off'); %direction selective to gratings
    end
    p_dir = find(p_anova_dir<0.05);


    % Do all fits at once
    [DSIstruct, ZpZcStruct, plaid_corr, gratingFitStruct, ZpZcPWdist, phaseModStruct] = bigFits(avg_resp_dir);

    % Get direction selectivity
        DSI         = DSIstruct.DSI;
        DSI_ind     = DSIstruct.DS_ind;
        DSI_maxInd  = DSIstruct.prefDir;
        g_dsi       = DSIstruct.gDSI;
        g_osi       = DSIstruct.gOSI;
        ang_dir     = DSIstruct.gDSI_prefDir;
        ang_ori     = DSIstruct.gOSI_prefDir;
    
    % Get direction tuning curve fit
        dir_b_hat_all           = gratingFitStruct.b;
        k1_hat_all          = gratingFitStruct.k1;
        R1_hat_all          = gratingFitStruct.R1;
        R2_hat_all          = gratingFitStruct.R2;
        u1_hat_all          = gratingFitStruct.u1;
        u2_hat_all          = gratingFitStruct.u2;
        dir_sse_all         = gratingFitStruct.sse;
        dir_R_square_all    = gratingFitStruct.Rsq;
        dir_yfit_all        = gratingFitStruct.yfit;

    % Get partial correlations
        Zp = ZpZcStruct.Zp;
        Zc = ZpZcStruct.Zc;
        Rp = ZpZcStruct.Rp;
        Rc = ZpZcStruct.Rc;
    % Get PCI fit, get amplitude and baseline
        PCI             = phaseModStruct.PCI;
        yfit_all        = phaseModStruct.yfit;
        amp_hat_all     = phaseModStruct.amp;
        b_hat_all       = phaseModStruct.b;
        sse_all         = phaseModStruct.sse;
        R_square_all    = phaseModStruct.rsq;

% 
% 
% % global DSI
%     angs = 0:30:330;
%     eps_val = 1e-3;  % Small epsilon to prevent division by zero
%     for j = 1:nCells
%         amps        = avg_resp_dir(j,:,1,1,1); %our responses  
%         amps(amps < 0) = 0;
%         total_response = sum(amps); % Normalize factor
%         if total_response < eps_val
%             g_dsi(j) = NaN; % Assign NaN if total response is too small
%             g_osi(j) = NaN;
%             ang(j) = NaN;
%             ang_ori(j) = NaN;
%         continue
%         end
%         % Direction Selectivity Index (gDSI)
%             vec_x = sum(cos(deg2rad(angs)) .* amps);
%             vec_y = sum(sin(deg2rad(angs)) .* amps);
%             g_dsi(j) = sqrt(vec_x^2 + vec_y^2) / total_response;
%         % Orientation Selectivity Index (gOSI)
%             vec_x_ori = sum(cos(deg2rad(2 * angs)) .* amps);
%             vec_y_ori = sum(sin(deg2rad(2 * angs)) .* amps);
%             g_osi(j) = sqrt(vec_x_ori^2 + vec_y_ori^2) / total_response;
%         % Preferred direction (in degrees)
%             ang_dir(j) = mod(rad2deg(atan2(vec_y, vec_x)), 360);
%         % Preferred orientation (in degrees, folded into [0, 180))
%             ang_ori(j) = mod(0.5 * rad2deg(atan2(vec_y_ori, vec_x_ori)), 180);
%     end
% 

    if ~exist(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)]))
        mkdir(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)]))
    end
    save(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_fitsSG.mat']), 'g_dsi', 'g_osi', 'ang_dir', 'ang_ori', 'resp_ind_dir','nCells','nTrials','nDirs','avg_resp_dir','p_dir','DSI','plaid_corr','Rp','Rc','Zp','Zc','ZpZcPWdist','yfit_all','amp_hat_all','b_hat_all','sse_all','R_square_all','dir_yfit_all','k1_hat_all','dir_sse_all','dir_R_square_all');


%% set inclusion criteria
resp_ind    = intersect(intersect(resp_ind_dir,p_dir),find(DSI>0.5));
ind         = resp_ind;


if doPlot == 1
%% Plot population ZpZc
plotZpZc4PhasePopulation(ZpZcStruct,resp_ind,30)
    movegui('center')
    sgtitle('Pattern direction selective cells at four phases')
print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_ZpZc.pdf']),'-dpdf', '-fillpage') 


%% Plot grating tuning curves by cell


ind=resp_ind;

if doGratPlot == 1

respgrat = avg_resp_dir(:,:,1,1,1);
semgrat = avg_resp_dir(:,:,1,1,2);
[prefresp, prefdir] = max(respgrat,[],2);
vecshift = zeros(size(respgrat,1),1)+6; 
vecshift = vecshift - prefdir;

avg_resp_grat =[];
std_resp_grat = [];


for iCell = 1:size(respgrat,1)
    avg_resp_grat(iCell,:) = circshift(respgrat(iCell,:),vecshift(iCell),2);
    sem_resp_grat(iCell,:) = circshift(semgrat(iCell,:),vecshift(iCell),2);
end

figure;
start = 1;
n = 1;
x=[-150:30:180];
for iCell = 1:nCells
    subplot(5,4,start)
        plot(x, avg_resp_grat(iCell,:),'k') 
        shadedErrorBar(x,avg_resp_grat(iCell,:),sem_resp_grat(iCell,:))
        xlabel('grating direction')
        ylabel('sp/s')
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell)]); 
        end
    start = start+1;
    if start >20
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Direction tuning curves aligned to preferred grating direction'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_DirectionTuningGratings_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Direction tuning curves aligned to preferred grating direction'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_DirectionTuningGratings_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end
end
close all

else
end

%% Plot pattern index modulation by cell

PCI = (Zp-Zc);

%fit sinusoid
phase = [0 90 180 270];
phase_range = 0:1:359;
figure;
start=1;
n=1;
for iCell = 1:length(ind)
    ic=ind(iCell);
    subplot(5,4,start)
        scatter(phase,PCI(:,ic),8,'filled');
        hold on
        [b_hat_all(ic,1), amp_hat_all(ic,1), per_hat_all(ic,1),pha_hat_all(ic,1),sse_all(ic,1),R_square_all(ic,1)] = sinefit_PCI(deg2rad(phase),PCI(:,ic));
        yfit_all(ic,:,1) = b_hat_all(ic,1)+amp_hat_all(ic,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(ic,1) + 2.*pi/pha_hat_all(ic,1)));
        idx = ic==ind; %does ic equal any value in the index?
        if any(idx)  %if yes, plot in red. if not, plot in black
            plot(phase_range, yfit_all(ic,:,1),'k'); 
            subtitle(['cell ' num2str(ic) ', Rsq ' num2str(R_square_all(ic),'%.3f'), ', SSE ' num2str(sse_all(ic),'%.2f')],'fontweight','bold')
        else
            plot(phase_range, yfit_all(ic,:,1),'k:');
            subtitle(['cell ' num2str(ic) ', Rsq ' num2str(R_square_all(ic),'%.2f'), ', SSE ' num2str(sse_all(ic),'%.2f')])
        end
        ylabel('Zp-Zc'); xlabel('Mask phase'); ylim([-7 7])
        xlim([0 360]); xticks([0 180 360]); set(gca,'TickDir','out'); axis square
    start = start+1;
    if start >20
        sgtitle(['marmosetV1_' expts(iexp,:) ' PCI modulation across mask phase by cell'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_PCImodulation_' num2str(n) '.pdf']),'-dpdf', '-fillpage')       
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle(['marmosetV1_' expts(iexp,:) ' PCI modulation across mask phase by cell'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_PCImodulation_' num2str(n) '.pdf']),'-dpdf', '-fillpage')  
    end        
end
close all


%% Plot polar plots by cell

figure;
start = 1;
n = 1;

[avg_resp_grat avg_resp_plddir] = getAlignedGratPlaidTuning(avg_resp_dir);

x=[-150:30:180];
x_rad = deg2rad(x);
for iCell =1:length(ind)
    ic = ind(iCell);
    subplot(5,4,start)
        for im = 1:nPhas
            polarplot([x_rad x_rad(1)], [avg_resp_plddir(ic,:,im) avg_resp_plddir(ic,1,im)])
            hold on
        end
        polarplot([x_rad x_rad(1)], [avg_resp_grat(ic,:) avg_resp_grat(ic,1)],'k', 'LineWidth',2) 
        idx = ic==ind;
        if any(idx)
            subtitle(['cell ' num2str(ic)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(ic)])
        end
    start = start+1;    
    if start>20
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Polar plots'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_PolarPlots_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Polar plots'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_PolarPlots_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end
end     
close all

%% Plot Zp Zc classifications by cell


sz  = 25;

figure;
start = 1;
n = 1;
for iCell = 1:nCells
    subplot(5,4,start)
        plotZpZc4PhaseCell(ZpZcStruct,iCell,sz)
        if iCell ==1; legend('0 deg','90 deg','180 deg', '270 deg'); end;
        idx = iCell==ind;
        if any(idx)
            subtitle(['cell ' num2str(iCell)],'fontweight','bold')
        else
            subtitle(['cell ' num2str(iCell)])
        end
        plotZcZpBorders; set(gca,'TickDir','out'); axis square    
    start = start +1;
    if start>20
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Zp Zc by cell'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_ZcZpByCell_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        figure;
        movegui('center')
        start = 1;
        n = n+1;
    end
    if iCell == nCells
        sgtitle(['expt ' num2str(expts(iexp,:)) ' - Zp Zc by cell'])
        print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_ZcZpByCell_' num2str(n) '.pdf']), '-dpdf','-fillpage')
    end        
end
close all



%% Direction tuning curves with sem

[sem_resp_grat] = getAlignedGratSEM(avg_resp_dir);

gDSIvals = round(g_dsi,2);


figure;
n=1;
start=1;
x=[-150:30:180];
    for iCell = 1:nCells
        subplot(6,3,start)
            plot(x, avg_resp_grat(iCell,:),'k') 
            hold on
            shadedErrorBar(x,avg_resp_grat(iCell,:),sem_resp_grat(iCell,:))
            ylim([0 20])
            ylabel('hz')
            DirIdx = iCell==p_dir;
            IndIdx = iCell==ind;
            if any(IndIdx)
                subtitle([num2str(iCell) ' - gDSI=' num2str(gDSIvals(iCell)) ' anova ✓'],'fontweight','bold')
            elseif any(DirIdx)
                subtitle([num2str(iCell) ' - gDSI=' num2str(gDSIvals(iCell)) ' anova ✓'])
            else
                subtitle([num2str(iCell) ' - gDSI=' num2str(gDSIvals(iCell))])
            end
        start = start+1;    
        if start>18
            movegui('center')
            sgtitle(['expt ' num2str(expts(iexp,:)) ' - Grating tuning'])
            print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_GratingTuningCurve_' num2str(n) '.pdf']), '-dpdf','-fillpage')
            figure;
            movegui('center')
            start = 1;
            n = n+1;
        end
        if iCell == nCells
            movegui('center')
            sgtitle(['expt ' num2str(expts(iexp,:)) ' - Grating tuning'])
            print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)], [svName '_marmosetV1_' expts(iexp,:) '_GratingTuningCurve_' num2str(n) '.pdf']), '-dpdf','-fillpage')
        end
    end
   




else 
end


