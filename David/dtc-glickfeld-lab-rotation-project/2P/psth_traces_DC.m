%% Per-Cell PSTH — raw individual trial df/f traces, color-coded by condition
% Shows each trial's trace (baseline + stimulus period) aligned to stimulus
% onset. Grating=black, plaid=blue, blank=gray. Paginated 5x4 (20 cells/page).

clc; clear all; close all;
startup
ds = 'CrossOriRandDir_ExptList_DC';
eval(ds)
iexp = 3;
frame_rate = 15;

%% Path construction & data loading
mouse = expt(iexp).mouse;
date = expt(iexp).date;
ImgFolder = expt(iexp).coFolder;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

dataBase = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\' expt(iexp).saveLoc];
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\David';

sessionDir = fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);

fprintf([mouse ' ' date '\n'])
fprintf('Loading data...\n')
load(fullfile(sessionDir, [date '_' mouse '_' run_str '_TCs.mat']))
load(fullfile(sessionDir, [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(sessionDir, [date '_' mouse '_' run_str '_input.mat']))

%% Compute trial-aligned df/f
prewin_frames = unique(celleqel2mat_padded(input.tItiWaitFrames))./3;
postwin_frames = unique(celleqel2mat_padded(input.nStimOneFramesOn));
tt = (1-prewin_frames:postwin_frames).*(1/frame_rate);

data_resp = nan(prewin_frames+postwin_frames, nCells, nTrials);
for itrial = 1:nTrials
    if cStimOn(itrial) + postwin_frames < sz(3)
        data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
    end
end

data_f = mean(data_resp(1:prewin_frames,:,:), 1);
data_dfof = (data_resp - data_f) ./ data_f;

%% Condition indices
ind_stimAlone = intersect(find(stimCon_all), find(maskCon_all==0));
ind_plaid     = intersect(find(stimCon_all), find(maskCon_all));
ind_blank     = intersect(find(stimCon_all==0), find(maskCon_all==0));

%% Paginated per-cell PSTH plot (20 cells/page, 5x4 grid)
outDir = fullfile(base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], 'dfof_traces');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

figure('Visible','off', 'Position', [0 0 1920 1080]);
start = 1;
n = 1;
for iCell = 1:nCells
    subplot(5,4,start)
    hold on

    % Blank trials (gray, behind)
    for it = 1:length(ind_blank)
        plot(tt, data_dfof(:,iCell,ind_blank(it)), 'Color', [0.6 0.6 0.6 0.2], 'LineWidth', 0.3)
    end
    % Grating trials (black)
    for it = 1:length(ind_stimAlone)
        plot(tt, data_dfof(:,iCell,ind_stimAlone(it)), 'Color', [0 0 0 0.2], 'LineWidth', 0.3)
    end
    % Plaid trials (blue)
    for it = 1:length(ind_plaid)
        plot(tt, data_dfof(:,iCell,ind_plaid(it)), 'Color', [0 0 1 0.2], 'LineWidth', 0.3)
    end

    xline(0, 'r--')
    subtitle(['cell ' num2str(iCell)])
    if start == 1
        xlabel('time (s)')
        ylabel('df/f')
    end
    set(gca, 'FontSize', 5)
    hold off

    start = start + 1;
    if start > 20
        sgtitle({[mouse ' ' date ' — Per-cell PSTH: raw trial traces'], ...
            ['grating (black) ' num2str(length(ind_stimAlone)) 't, ' ...
             'plaid (blue) ' num2str(length(ind_plaid)) 't, ' ...
             'blank (gray) ' num2str(length(ind_blank)) 't']}, 'FontSize', 10)
        print(fullfile(outDir, [date '_' mouse '_' run_str '_PSTH_' num2str(n) '.png']), '-dpng', '-r300')
        figure('Visible','off', 'Position', [0 0 1920 1080]);
        start = 1;
        n = n + 1;
    end
    if iCell == nCells
        sgtitle({[mouse ' ' date ' — Per-cell PSTH: raw trial traces'], ...
            ['grating (black) ' num2str(length(ind_stimAlone)) 't, ' ...
             'plaid (blue) ' num2str(length(ind_plaid)) 't, ' ...
             'blank (gray) ' num2str(length(ind_blank)) 't']}, 'FontSize', 10)
        print(fullfile(outDir, [date '_' mouse '_' run_str '_PSTH_' num2str(n) '.png']), '-dpng', '-r300')
    end
end

fprintf('Saved %d pages to %s\n', n, outDir)
