%% Analyze raw marmoset data sent from Nicholas Priebe, 4/18/2025
clc; clear all; close all;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
svName = 'randDirFourPhase_CrossOri';
bnsz = 10; %spike binsize = 10 ms

doPlot                  = 1;
doGratPlot              = 1;
doBaselineSubtraction   = 1;
doPeakResponse          = 1;

expts = strvcat('g01', 'g06b', 'g12b', 'g17b','tss2','tss4','tss6','tss7','elf1','elf3'); %Four acute penetrations in V1, 1 marmoset
nexp = length(expts);
% iexp = 1;
% expt = expts(iexp,:);

%'g01', 'g06', 'g12', 'g17','tss2','tss4','tss6','tss7','elf1','elf3'
expt = 'g12';

    fprintf([expt '\n'])

    % Load data from experiment. Each experiment file contains one
    % variable, resp, that is nCells by nDir by nPhase by Grating/Plaid by
    % nTrials by Time (in 10 ms bins). For example, 310x12x4x2x13x120.
    % Grating type only has the first phase filled in.
    load(fullfile(base, 'Data\fromNicholas\CrossOri_randDirFourPhase_V1_marmoset', ([expt '.mat'])))
    
    nCells = size(resp,1);
    nDirs = size(resp,2);
    nPhas = size(resp,3);
    nTrials = size(resp,5);
    nStim = (size(resp,3)+1)*size(resp,2);


    % resp_cell = resp(:,:,:,:,:,6:105);
    % base_cell = resp(:,:,:,:,:,1:5);
    resp_cell = resp(:,:,:,:,:,21:end);
    base_cell = resp(:,:,:,:,:,1:20);

    if doPeakResponse
        nT2 = 2*nTrials;
        nBinsCyc = 50;   % 500 ms per cycle at 10 ms bins
        binWidth = 0.01;   % seconds per bin

        % Split into 2 cycles and stack along the trial dimension (dim 5)
        % -> nCells x nDirs x nPhas x 2 x (2*nTrials) x 50
        cyc1 = resp_cell(:,:,:,:,:,1:nBinsCyc);
        cyc2 = resp_cell(:,:,:,:,:,nBinsCyc+1:2*nBinsCyc);
        resp_cyc = cat(5, cyc1, cyc2);
        
        % Outputs
        peak     = nan(nCells, nDirs, nPhas, 2);   % mean dc + f1
        f1_all   = nan(nCells, nDirs, nPhas, 2);
        dc_all   = nan(nCells, nDirs, nPhas, 2);
        peak_trl = nan(nCells, nDirs, nPhas, 2, nT2);  % per-trial peak (for error/stats)
        
        for iC = 1:nCells
          for iD = 1:nDirs
            for iP = 1:nPhas
              for iS = 1:2   % 1 = grating, 2 = plaid
                % nT2 x 50 matrix of single-cycle data
                data = squeeze(resp_cyc(iC,iD,iP,iS,:,:));
        
               % NaN-drop and peak_trl assignment
                valid = ~any(isnan(data),2);
                if ~any(valid), continue; end
                idx  = find(valid);            % positions in the 2*nTrials stack
                data = data(valid,:) / binWidth;
                
                [f1, f1s, dc, dcs] = compf1werr(data);
                
                f1_all(iC,iD,iP,iS) = f1;
                dc_all(iC,iD,iP,iS) = real(dc);
                peak(iC,iD,iP,iS)   = real(dc) + f1;
                peak_trl(iC,iD,iP,iS,idx) = reshape(real(dcs(:)) + f1s(:), [1 1 1 1 numel(idx)]);
                peak_trl_f1(iC,iD,iP,iS,idx) = reshape(f1s(:), [1 1 1 1 numel(idx)]);
              end
            end
          end
        end

        avg_resp_dir_f1dc(:,:,:,:,1) = mean(peak_trl,5, 'omitnan'); % Average across time and trial
        avg_resp_dir_f1dc(:,:,:,:,2) = std(peak_trl,0,5, 'omitnan')/sqrt(nT2); % SEM for each stim condition
        avg_resp_dir_f1(:,:,:,:,1) = mean(peak_trl_f1,5, 'omitnan'); % Average across time and trial
        avg_resp_dir = avg_resp_dir_f1dc;
    else

        avg_resp_dir(:,:,:,:,1) = mean(sum(resp_cell,6),5); % Average across time and trial
        avg_resp_dir(:,:,:,:,2) = std(sum(resp_cell,6),0,5)/sqrt(nTrials); % SEM for each stim condition
    end

    resp_dir_tc = mean(resp,5); % Average across trial, to look at time 
    resp_dir_tr = mean(resp,6); % Average across time, to look at trial 

    if doPeakResponse
        % per-trial peak (dc + f1, Hz), averaged over the 2 cycles of each trial
        resp_cell_trials = mean(cat(6, peak_trl(:,:,:,:,1:nTrials), peak_trl(:,:,:,:,nTrials+1:end)), 6, 'omitnan');
        % -> nCells x nDirs x nPhas x 2 x nTrials, same shape as before
    else
        resp_cell_trials = sum(resp_cell,6) * 1;   % mean Hz per trial
    end

    % Baseline stays a mean rate. Spontaneous activity has no 2 Hz modulation, and 0.2 s is less than one cycle, so an f1 fit isn't meaningful
    base_cell_trials = sum(base_cell,6) * 5; % Hz per trial

    mean_base_all = mean(base_cell_trials, [2 3 4 5], 'omitnan'); % mean across all stim conditions
    if doBaselineSubtraction
        avg_resp_dir(:,:,:,:,1) = avg_resp_dir(:,:,:,:,1) - mean_base_all;
    end


    
    for id = 1:nDirs
         resp_tmp = squeeze(resp_cell_trials(:,id,1,1,:));  % [nCells x nTrials]
         base_tmp = squeeze(base_cell_trials(:,id,1,1,:));
        [h_resp(:,id,1,1), p_resp(:,id,1,1)] = ttest(resp_tmp,base_tmp,'dim', 2,'tail', 'right', 'alpha', (0.05 / nDirs));
    end
    resp_ind_dir = find(sum(h_resp(:,:,1,1),2)); %sig responsive to gratings
    

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
    [DSIstruct, ZpZcStruct, plaid_corr, gratingFitStruct, ZpZcPWdist, phaseModStruct] = bigFits(avg_resp_dir,1); % 2nd input for marmoset data 

    DSIstruct = getDSIstruct_new(avg_resp_dir);
    % Get direction selectivity
        DSI         = DSIstruct.DSI;
        DSI_ind     = DSIstruct.DS_ind;
        DSI_maxInd  = DSIstruct.prefDir;
        g_dsi       = DSIstruct.gDSI;
        g_osi       = DSIstruct.gOSI;
        ang_dir     = DSIstruct.gDSI_prefDir;
        ang_ori     = DSIstruct.gOSI_prefDir;
    
    % Get direction tuning curve fit
        dir_b_hat_all       = gratingFitStruct.b;
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

    if ~exist(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expt]))
        mkdir(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expt]))
    end
    save(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expt], [svName '_marmosetV1_' expt '_fitsSG.mat']),'avg_resp_dir_f1dc', 'g_dsi', 'g_osi', 'ang_dir', 'ang_ori', 'resp_ind_dir','nCells','nTrials','nDirs','avg_resp_dir','p_dir','DSI','plaid_corr','Rp','Rc','Zp','Zc','ZpZcPWdist','yfit_all','amp_hat_all','b_hat_all','sse_all','R_square_all','dir_yfit_all','k1_hat_all','dir_sse_all','dir_R_square_all');


%% set inclusion criteria
resp_ind    = intersect(resp_ind_dir,find(DSI>0.5));

max_grat_resp = max(avg_resp_dir(:,:,1,1,1), [], 2);
grat_ind = find(max_grat_resp>1);

ind = intersect(grat_ind, resp_ind_dir);

if doPlot == 1
%% Plot population ZpZc
plotZpZc4PhasePopulation(ZpZcStruct,ind,30)
    movegui('center')
    sgtitle('Pattern direction selective cells at four phases')
print(fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expt], [svName '_marmosetV1_' expt '_ZpZc.pdf']),'-dpdf', '-fillpage') 


%% Plot grating tuning curves by cell


% ind=resp_ind;

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

start = 1;
n = 1;

[avg_resp_grat, avg_resp_plddir] = getAlignedGratPlaidTuning(avg_resp_dir,'whole-cell');

% ind_Nicholas = [19, 27, 30, 37, 40, 41, 56, 67, 83, 84, 88, 89, 92, 93,95, 98, 99, 100, 102, 103,107, 111, 114, 121, 122, 124, 128, 129, 130, 134, 135, 137, 140, 141, 143, 147, 148, 149, 153, 154, 155, 158, 159, 160, 163, 165, 169, 170, 171, 172, 174, 175, 176, 177, 179, 181, 182, 187, 188, 189, 192, 202, 213, 214, 215, 220, 224, 227];

pdfDir = fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)]);
if ~exist(pdfDir, 'dir'); mkdir(pdfDir); end

pdfFile = fullfile(pdfDir, ['randDirFourPhase_CrossOri_marmosetV1_' expts(iexp,:) '_PolarPlots.pdf']);
if isfile(pdfFile); delete(pdfFile); end

pageSize = [0 0 8.5 11];      % letter, inches (use [0 0 11 8.5] for landscape)
nRows = 5; nCols = 4; perPage = nRows*nCols;

x = -150:30:180;
x_rad = deg2rad(x);
start = 1;
fig = [];

for iCell = 1:nCells
    ic = iCell;

    % New page: create a figure that fills the paper
    if start == 1
        fig = figure('Units','inches','Position',pageSize, ...
                     'PaperUnits','inches','PaperSize',pageSize(3:4), ...
                     'PaperPosition',pageSize,'Color','w');
    end

    % Tight subplot layout (subplot's own spacing wastes a lot of room)
    row = ceil(start/nCols);
    col = mod(start-1, nCols) + 1;
    mL = 0.03; mR = 0.03; mB = 0.03; mT = 0.05;   % page margins (normalized)
    gap = 0.01;
    w = (1 - mL - mR - gap*(nCols-1)) / nCols;
    h = (1 - mT - mB - gap*(nRows-1)) / nRows;
    scale = 0.6;   % 0.8 = 20% smaller; adjust to taste
    % Cell (slot) position on the page
    cellPos = [mL + (col-1)*(w+gap), 1 - mT - row*h - (row-1)*gap, w, h];
    % Shrink the axes about the slot center
    pos = [cellPos(1) + cellPos(3)*(1-scale)/2, ...
           cellPos(2) + cellPos(4)*(1-scale)/2, ...
           cellPos(3)*scale, cellPos(4)*scale];

    pax = polaraxes('Parent',fig,'Position',pos);
    pax.TitleFontWeight = 'normal';

    hold(pax,'on')
    for im = 1:nPhas
        polarplot(pax,[x_rad x_rad(1)], [avg_resp_plddir(ic,:,im) avg_resp_plddir(ic,1,im)])
    end
    polarplot(pax,[x_rad x_rad(1)], [avg_resp_grat(ic,:) avg_resp_grat(ic,1)],'k','LineWidth',2)

    isSel = any(ic == ind);
    if isSel
        title(pax, ['cell ' num2str(ic)], 'FontWeight','bold','FontSize',9)
       
    else
        title(pax, ['cell ' num2str(ic)], 'FontSize',9)
    end
    % isSelN = any(ic == ind_Nicholas);
    % if isSelN
    %      % Black box around the whole subplot cell
    %     annotation(fig,'rectangle',[pos(1)-gap/2, pos(2)-gap/2, pos(3)+gap, pos(4)+gap], ...
    %                'Color','k','LineWidth',1.5)
    % end

    start = start + 1;

    % Page full, or last cell: export and close
    if start > perPage || iCell == nCells
        exportgraphics(fig, pdfFile, 'ContentType','vector','Append',true);
        close(fig)
        start = 1;
    end
end
close all



figure;
    subplot 221
            histogram(base_cell_trials(30,:,1,1,:),30); hold on
            histogram(resp_cell_trials(30,:,1,1,:),30)
            subtitle('cell 30')
    subplot 222
            histogram(base_cell_trials(37,:,1,1,:),30); hold on
            histogram(resp_cell_trials(37,:,1,1,:),30)
            subtitle('cell 37')
    subplot 223
            histogram(base_cell_trials(34,:,1,1,:),30); hold on
            histogram(resp_cell_trials(34,:,1,1,:),30)
            subtitle('cell 34')
    subplot 224
            histogram(base_cell_trials(40,:,1,1,:),30); hold on
            histogram(resp_cell_trials(40,:,1,1,:),30)
            subtitle('cell 40')

%% Plot tuning curves by cell (Cartesian)

start = 1;
n = 1;

[avg_resp_grat, avg_resp_plddir] = getAlignedGratPlaidTuning(avg_resp_dir,'whole-cell');

% ind_Nicholas = [19, 27, 30, 37, 40, 41, 56, 67, 83, 84, 88, 89, 92, 93,95, 98, 99, 100, 102, 103,107, 111, 114, 121, 122, 124, 128, 129, 130, 134, 135, 137, 140, 141, 143, 147, 148, 149, 153, 154, 155, 158, 159, 160, 163, 165, 169, 170, 171, 172, 174, 175, 176, 177, 179, 181, 182, 187, 188, 189, 192, 202, 213, 214, 215, 220, 224, 227];

pdfDir = fullfile(base, 'Analysis\Neuropixel\marmosetFromNicholas', ['marmosetV1_' expts(iexp,:)]);
if ~exist(pdfDir, 'dir'); mkdir(pdfDir); end

pdfFile = fullfile(pdfDir, ['randDirFourPhase_CrossOri_marmosetV1_' expts(iexp,:) '_TuningPlots.pdf']);
if isfile(pdfFile); delete(pdfFile); end

pageSize = [0 0 8.5 11];      % letter, inches (use [0 0 11 8.5] for landscape)
nRows = 5; nCols = 4; perPage = nRows*nCols;

x = -150:30:180;              % direction (deg)
start = 1;
fig = [];

for iCell = 1:nCells
    ic = iCell;

    % New page: create a figure that fills the paper
    if start == 1
        fig = figure('Units','inches','Position',pageSize, ...
                     'PaperUnits','inches','PaperSize',pageSize(3:4), ...
                     'PaperPosition',pageSize,'Color','w');
    end

    % Tight subplot layout
    row = ceil(start/nCols);
    col = mod(start-1, nCols) + 1;
    mL = 0.05; mR = 0.03; mB = 0.05; mT = 0.05;   % page margins (normalized)
    gap = 0.03;                                    % a bit larger for tick labels
    w = (1 - mL - mR - gap*(nCols-1)) / nCols;
    h = (1 - mT - mB - gap*(nRows-1)) / nRows;
    scale = 0.7;   % shrink axes about slot center; adjust to taste
    cellPos = [mL + (col-1)*(w+gap), 1 - mT - row*h - (row-1)*gap, w, h];
    pos = [cellPos(1) + cellPos(3)*(1-scale)/2, ...
           cellPos(2) + cellPos(4)*(1-scale)/2, ...
           cellPos(3)*scale, cellPos(4)*scale];

    ax = axes('Parent',fig,'Position',pos);
    ax.TitleFontWeight = 'normal';
    ax.FontSize = 7;
    ax.Box = 'off';

    hold(ax,'on')
    for im = 1:nPhas
        plot(ax, x, squeeze(avg_resp_plddir(ic,:,im)))
    end
    plot(ax, x, avg_resp_grat(ic,:), 'k', 'LineWidth', 2)

    xlim(ax, [x(1) x(end)])
    xticks(ax, -150:90:180)
    if row == nRows
        xlabel(ax, 'Direction (deg)', 'FontSize', 7)
    end
    if col == 1
        ylabel(ax, 'Response', 'FontSize', 7)
    end

    isSel = any(ic == ind);
    if isSel
        title(ax, ['cell ' num2str(ic)], 'FontWeight','bold','FontSize',9)
    else
        title(ax, ['cell ' num2str(ic)], 'FontSize',9)
    end

    % isSelN = any(ic == ind_Nicholas);
    % if isSelN
    %     % Black box around the whole subplot cell
    %     annotation(fig,'rectangle', ...
    %                [cellPos(1)+gap/4, cellPos(2)-gap/4, cellPos(3)-gap/2, cellPos(4)+gap/4], ...
    %                'Color','k','LineWidth',1.5)
    % end

    start = start + 1;

    % Page full, or last cell: export and close
    if start > perPage || iCell == nCells
        exportgraphics(fig, pdfFile, 'ContentType','vector','Append',true);
        close(fig)
        start = 1;
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
                subtitle([num2str(iCell) ' - gDSI=' num2str(gDSIvals(iCell)) ' anova '],'fontweight','bold')
            elseif any(DirIdx)
                subtitle([num2str(iCell) ' - gDSI=' num2str(gDSIvals(iCell)) ' anova '])
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


   
%% plot not included cells
%%%%%%%%%%%%%%%%%%%%%%%%

ind_DS = intersect(find(DSI>0.25),find(g_dsi>0.2));
    [DSIstruct_new] = getDSIstruct_new(avg_resp_dir);
    gratResp        = DSIstruct_new.resp;
    maxGratResp     = max(gratResp,[],2);
    minGratResp     = min(gratResp,[],2);
ind_peakFRmin   = find(maxGratResp>1);
ind_rangeFRmin  = find(maxGratResp > (abs(minGratResp)*2));
% ind_FR = intersect(ind_peakFRmin,ind_rangeFRmin);
ind_FR = ind_peakFRmin;



% get pref dir
gratResps           = squeeze(avg_resp_dir(:,:,1,1,1));
[maxDir, indDir]    = max(gratResps,[],2);
dirs                = 0:30:330;

% get aligned tuning curves
[avg_resp_grat, avg_resp_plaid, sem_resp_grat, sem_resp_plaid] = getAlignedGratPlaidTuning(avg_resp_dir);
[avg_resp_grat_f1dc, avg_resp_plaid_f1dc, ~, ~] = getAlignedGratPlaidTuning(avg_resp_dir_f1dc,'','avgonly');

x       = [-150:30:180];
x_rad   = deg2rad(x);

% for PCI fit
PCI = (Zp-Zc);
phase = [0 90 180 270];
phase_range = 0:1:359;

% get colors
colors  = getColors;


% Print STA time point choices
pdfDir = fullfile(base, 'Analysis', 'Neuropixel','marmosetFromNicholas', ['marmosetV1_' expt]);

pdfFile = fullfile(pdfDir, [expt '-allCells_withTuning.pdf']);
if isfile(pdfFile); delete(pdfFile); end

for ic = 1:nCells
    iCell = ic;

    figure();
    sgtitle(sprintf('cell %d', iCell), 'FontSize', 10); 

    subplot(3,5,1)
        plot(x, avg_resp_grat_f1dc(iCell,:),'k')
        subtitle('grating tuning','FontSize', 6)
        xlabel('direction','FontSize', 6)
        ylabel('Hz (f1+dc)','FontSize', 6)
        subtitle('F1+DC')
        set(gca,'FontSize',5,'TickDir','out'); box off
    subplot(3,5,2)
        plot(x, avg_resp_grat(iCell,:),'k')
        subtitle('grating tuning','FontSize', 6)
        xlabel('direction','FontSize', 6)
        ylabel('Hz (f1+dc, bl sub)','FontSize', 6)
        subtitle('F1+DC, baseline subtracted')
        set(gca,'FontSize',5,'TickDir','out'); box off
    subplot(3,5,3)
        data = avg_resp_grat(iCell,:);
        [minVal, ~] = min(data);
            if minVal < 0
                resp = data-minVal;
            else
                resp = data;
            end
        plot(x, resp, 'k')
        subtitle(['gDSI=' num2str(round(g_dsi(iCell),2)) ', DSI=' num2str(round(DSI(iCell),2))],'FontSize', 6)
        xlabel('direction','FontSize', 6)
        ylabel('Hz (used only for DSI)','FontSize', 6)
        set(gca,'FontSize',5,'TickDir','out'); box off
    subplot(3,5,4)
        for im = 1:4
            polarplot([x_rad x_rad(1)], [avg_resp_plaid_f1dc(iCell,:,im) avg_resp_plaid_f1dc(iCell,1,im)],'Color', colors(im,:))
            hold on
        end
        polarplot([x_rad x_rad(1)], [avg_resp_grat_f1dc(iCell,:) avg_resp_grat_f1dc(iCell,1)],'k', 'LineWidth',2) 
        set(gca,'FontSize',5)
    subplot(3,5,5)
        for im = 1:4
            scatter(Zc(im,iCell), Zp(im,iCell),8,colors(im,:),'filled')
            hold on
        end
        ylabel('Zp'); ylim([-4 8]);
        xlabel('Zc'); xlim([-4 8]);
        plotZcZpBorders; axis square; set(gca,'TickDir','out')
    subplot(3,5,10)
        scatter(phase,PCI(:,iCell),8,'filled'); hold on
        [b_hat_all(iCell,1), amp_hat_all(iCell,1), per_hat_all(iCell,1),pha_hat_all(iCell,1),sse_all(iCell,1),R_square_all(iCell,1)] = sinefit_PCI(deg2rad(phase),PCI(:,iCell));
        yfit_all(iCell,:,1) = b_hat_all(iCell,1)+amp_hat_all(iCell,1).*(sin(2*pi*deg2rad(phase_range)./per_hat_all(iCell,1) + 2.*pi/pha_hat_all(iCell,1)));
        plot(phase_range, yfit_all(iCell,:,1),'k:');
        ylabel('Zp-Zc'); xlabel('Mask phase'); ylim([-7 7])
        xlim([0 360]); xticks([0 180 360]); set(gca,'TickDir','out')
        text(0.05,0.9,sprintf('R^2 = %.2f',R_square_all(iCell,1)),'Units','normalized','FontSize',5) 

    subplot(3,5,6)
        axis off
            txtOpts = {'Units','normalized','Color','k','FontSize',5,'HorizontalAlignment','left','VerticalAlignment','top'};
            txtOptsBold = {'Units','normalized','Color','k','FontSize',5,'FontWeight','bold','HorizontalAlignment','left','VerticalAlignment','top'};
            if ismember(ic,resp_ind_dir)
                text(0.02, 0.98, 'pass Vis Resp criteria', txtOptsBold{:});
            else
                text(0.02, 0.98, 'pass Vis Resp criteria', txtOpts{:});
            end
            if ismember(ic,ind_DS)
                text(0.02, 0.68, 'pass DS criteria - gDSI>0.2, DSI>0.25', txtOptsBold{:});
            else
                text(0.02, 0.68, 'pass DS criteria - gDSI>0.2, DSI>0.25', txtOpts{:});
            end
            if ismember(ic,ind_FR)
                text(0.02, 0.38, 'pass FR criteria - peak grating FR > 1Hz', txtOptsBold{:});
            else
                text(0.02, 0.38, 'pass FR criteria - peak grating FR > 1Hz', txtOpts{:});
            end

    % Append current figure as a new page in the PDF
    exportgraphics(gcf, pdfFile,'ContentType', 'vector','Append', true);
    close(gcf)
end


end
