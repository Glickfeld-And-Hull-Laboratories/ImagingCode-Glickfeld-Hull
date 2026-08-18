%% aggregateSummary_gDSI_DC.m
% Loop over experiments, run per-session gDSI analysis pipeline, concatenate
% cells across sessions, save a pooled summary .mat, and produce population
% gDSI scatter + polar histogram figures.
%
% Replicates the core trial-alignment, condition-sorting, and gDSI steps from
% randDir_CrossOri_exptAnalysis_DC_gDSI.m for each session, then concatenates.
%
% Experiments are grouped by plaid_angle (from ExptList) so each angle gets
% its own output folder (aggregate_90/, aggregate_45/, etc.).
%
% Output per angle group: aggregate_XX/
%   aggregate_crossOriSummary_gDSI.mat  (saved with -struct)
%     S.stim_dirs, S.plaid_offset_deg
%     S.resp_mean  (nTotal x 2*nStimDir)  — [grat_dir1..N, plaid_dir1..N]
%     S.resp_sem   (nTotal x 2*nStimDir)
%     S.gDSI_grat, S.gDSI_plaid, S.prefDir_grat_deg, S.prefDir_plaid_deg
%     S.resp_to_grat, S.resp_to_any, S.anova_grat, S.anova_plaid, S.DS_grat, S.DS_plaid, S.DS_both
%     S.DSI_grat, S.dist_VS, S.dist_IOC
%     S.session_id, S.sessions (per-session metadata)

clc; clear all; close all; clear all global;
startup
ds = 'CrossOriRandDir_ExptList_DC';
eval(ds)

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\David';
addpath(genpath('Z:\All_Staff\home\David\repositories\dtc-glickfeld-lab-rotation-project\CrossOri_PatternV1_SG_DC'));

%% ====== USER SETTINGS — EDIT HERE ======
iexp_start   = 3;                          % range start
iexp_end     = length(expt);               % range end
iexp_include = [];                         % additional indices beyond range
iexp_exclude = [];                         % indices to remove
iexp_list = setdiff(union(iexp_start:iexp_end, iexp_include), iexp_exclude);
gDSI_threshold = 0.3;                     % DS cutoff for gDSI (both grat & plaid)
% =========================================

%% Group by plaid angle
all_angles = arrayfun(@(i) expt(i).plaid_angle, iexp_list);
unique_angles = unique(all_angles);
fprintf('Found %d plaid angle group(s): %s\n', length(unique_angles), mat2str(unique_angles))

% Storage for cross-angle KDE overlay figure
kde_delta_by_angle = cell(1, length(unique_angles));
kde_angle_labels   = zeros(1, length(unique_angles));
kde_n_DS           = zeros(1, length(unique_angles));

for ai = 1:length(unique_angles)
    this_angle = unique_angles(ai);
    angle_expts = iexp_list(all_angles == this_angle);

    if isempty(angle_expts)
        fprintf('\nSkipping plaid angle %d deg — no experiments after filtering.\n', this_angle)
        continue
    end

    fprintf('\n########## Plaid angle = %d deg (expts: %s) ##########\n', ...
        this_angle, mat2str(angle_expts))

    % --- Output directory ---
    outDir = fullfile(base, 'Analysis\2P\aggregate', sprintf('aggregate_%d', this_angle));
    if ~exist(outDir, 'dir'), mkdir(outDir); end

    %% Preallocate accumulators
    all_resp_mean    = [];
    all_resp_sem     = [];
    all_gDSI_grat    = [];
    all_gDSI_plaid   = [];
    all_prefDir_grat = [];
    all_prefDir_plaid = [];
    all_resp_flag    = [];
    all_anova_grat   = [];
    all_anova_plaid  = [];
    all_session_id   = [];
    all_Zp_VS        = [];
    all_Zp_IOC       = [];
    all_Zc_VSIOC     = [];
    all_DSI          = [];   % peak DSI (grating) per cell
    all_dist_VS      = [];   % angular distance to VS prediction
    all_dist_IOC     = [];   % angular distance to IOC prediction
    all_dist_C2      = [];   % angular distance to C2 prediction
    all_resp_any     = [];   % responsive to any direction (grating or plaid)
    any_nondegen     = false;  % track if any session has non-degenerate VS/IOC
    sessions = struct([]);
    cell_offset = 0;
    ref_stimDirs = [];  % to validate consistent directions across sessions

    %% Loop over experiments
    for idx = 1:length(angle_expts)
        iexp = angle_expts(idx);

        % --- Session info (same as exptAnalysis lines 12-21) ---
        mouse = expt(iexp).mouse;
        date  = expt(iexp).date;
        ImgFolder = expt(iexp).coFolder;
        time  = expt(iexp).coTime;
        nrun  = length(ImgFolder);
        run_str  = catRunName(cell2mat(ImgFolder), nrun);
        dataBase = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\' expt(iexp).saveLoc];

        fprintf('\n=== Session %d/%d: %s %s (expt %d) ===\n', idx, length(angle_expts), mouse, date, iexp)

        % --- Load crossOriScript outputs (same as exptAnalysis lines 30-32) ---
        fprintf('  Loading data...\n')
        sessDir = fullfile(dataBase, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);
        load(fullfile(sessDir, [date '_' mouse '_' run_str '_TCs.mat']))
        load(fullfile(sessDir, [date '_' mouse '_' run_str '_dataStim.mat']))
        load(fullfile(sessDir, [date '_' mouse '_' run_str '_input.mat']))

        fprintf('  nCells = %d, nTrials = %d, nStimDir = %d\n', nCells, nTrials, nStimDir)

        % --- Validate consistent direction grid across sessions ---
        if isempty(ref_stimDirs)
            ref_stimDirs = stimDirs;
        else
            if ~isequal(stimDirs, ref_stimDirs)
                error('Session %d (%s %s) has different stimDirs than session 1. Cannot concatenate.', iexp, mouse, date)
            end
        end

        % --- Derive experiment params (same as exptAnalysis lines 35-58) ---
        signed_offset = unique(maskDiffs(maskDiffs ~= 0));
        if isempty(signed_offset)
            signed_offset = 45;
        end
        plaid_offset = abs(signed_offset(1));
        offset_rad   = deg2rad(signed_offset(1));

        SF_val = expt(iexp).SF;
        if isfield(expt, 'TF_stim') && ~isempty(expt(iexp).TF_stim)
            speed_stim = expt(iexp).TF_stim / SF_val;
            speed_mask = expt(iexp).TF_mask / SF_val;
        else
            speed_stim = expt(iexp).TF / SF_val;
            speed_mask = speed_stim;
        end

        alpha_VS  = atan2(speed_mask * sin(offset_rad), speed_stim + speed_mask * cos(offset_rad));
        alpha_IOC = atan2(speed_mask - speed_stim * cos(offset_rad), speed_stim * sin(offset_rad));
        shift_VS_deg  = -rad2deg(alpha_VS);
        shift_IOC_deg = -rad2deg(alpha_IOC);

        % C2 (mask-tracking) prediction
        c2_angle_sess = signed_offset(1);       % C2 prediction on delta_pref polar plot
        shift_C2_deg_sess = -signed_offset(1);  % for scatter plots: y = x + shift_C2

        fprintf('  Plaid offset: %+.0f deg, speed_stim=%.1f, speed_mask=%.1f\n', signed_offset(1), speed_stim, speed_mask)
        fprintf('  VS shift: %.1f deg, IOC shift: %.1f deg, C2 deltapref: %+.1f deg\n', shift_VS_deg, shift_IOC_deg, c2_angle_sess)

        % --- Trial-aligned dF/F (same as exptAnalysis lines 66-80) ---
        fprintf('  Computing trial-aligned dF/F...\n')
        prewin_frames  = unique(celleqel2mat_padded(input.tItiWaitFrames)) ./ 3;
        nFramesOn      = unique(celleqel2mat_padded(input.nStimOneFramesOn));
        postwin_frames = nFramesOn;

        data_resp = nan(prewin_frames + postwin_frames, nCells, nTrials);
        for itrial = 1:nTrials
            if cStimOn(itrial) + postwin_frames < sz(3)
                data_resp(:,:,itrial) = npSub_tc(cStimOn(itrial)-prewin_frames : cStimOn(itrial)+postwin_frames-1, :);
            end
        end
        data_f       = mean(data_resp(1:prewin_frames,:,:), 1);
        data_dfof_tc = (data_resp - data_f) ./ data_f;

        % --- Guard: this script assumes single mask phase ---
        nMaskPhas = length(unique(maskPhas_all(maskCon_all > 0)));
        if nMaskPhas > 1
            error('Expt %d has %d mask phases. Aggregate script assumes 1 phase (pools all plaid trials per dir).', iexp, nMaskPhas);
        end

        % --- Sort trials by condition (same as exptAnalysis lines 82-137) ---
        fprintf('  Sorting trials into %d conditions...\n', nStimDir*2)
        ind_stimAlone = intersect(find(stimCon_all), find(maskCon_all == 0));
        ind_maskAlone = intersect(find(stimCon_all == 0), find(maskCon_all));
        ind_plaid     = intersect(find(stimCon_all), find(maskCon_all));
        ind_blank     = intersect(find(stimCon_all == 0), find(maskCon_all == 0));

        resp_win   = prewin_frames+5 : prewin_frames+nFramesOn;
        resp_blank = squeeze(mean(data_dfof_tc(resp_win,:,ind_blank), 1));  % nCells x nBlankTrials

        avg_resp    = zeros(nCells, nStimDir*2);   % [grat_1..N, plaid_1..N]
        avg_resp_se = zeros(nCells, nStimDir*2);
        h_resp_all  = zeros(nCells, nStimDir*2);

        all_resp_grat_trials  = [];   all_dir_labels   = [];
        all_resp_plaid_trials = [];   all_plaid_labels = [];

        for iDir = 1:nStimDir
            ind_stimdir = find(stimDir_all == stimDirs(iDir));
            ind_maskdir = find(maskDir_all == stimDirs(iDir));
            ind_diralone = [intersect(ind_stimdir, ind_stimAlone), intersect(ind_maskdir, ind_maskAlone)];
            ind_dirplaid = intersect(ind_stimdir, ind_plaid);

            % --- Grating condition (column iDir) ---
            resp_grat = squeeze(mean(data_dfof_tc(resp_win,:,ind_diralone), 1));  % nCells x nTrials
            avg_resp(:, iDir)    = mean(resp_grat, 2);
            avg_resp_se(:, iDir) = std(resp_grat, [], 2) ./ sqrt(size(resp_grat, 2));
            [h_resp_all(:,iDir), ~] = ttest2(resp_grat, resp_blank, 'dim', 2, 'tail', 'right', 'alpha', 0.05/nStimDir);

            all_resp_grat_trials = [all_resp_grat_trials, resp_grat];
            all_dir_labels       = [all_dir_labels, iDir*ones(1, size(resp_grat, 2))];

            % --- Plaid condition (column nStimDir + iDir) ---
            col = nStimDir + iDir;
            resp_plaid = squeeze(mean(data_dfof_tc(resp_win,:,ind_dirplaid), 1));  % nCells x nTrials (all phases pooled)
            avg_resp(:, col)    = mean(resp_plaid, 2);
            avg_resp_se(:, col) = std(resp_plaid, [], 2) ./ sqrt(size(resp_plaid, 2));
            [h_resp_all(:,col), ~] = ttest2(resp_plaid, resp_blank, 'dim', 2, 'tail', 'right', 'alpha', 0.05/nStimDir);

            all_resp_plaid_trials = [all_resp_plaid_trials, resp_plaid];
            all_plaid_labels      = [all_plaid_labels, iDir*ones(1, size(resp_plaid, 2))];
        end

        % --- Responsiveness flags (same as exptAnalysis lines 139-153) ---
        fprintf('  Computing responsiveness flags...\n')
        resp_flag = any(h_resp_all(:, 1:nStimDir), 2)';   % responsive to any grating dir
        resp_any_flag = any(h_resp_all, 2)';              % responsive to any dir (grat or plaid)

        anova_g = zeros(1, nCells);
        anova_p = zeros(1, nCells);
        for iCell = 1:nCells
            anova_g(iCell) = anova1(all_resp_grat_trials(iCell,:), all_dir_labels, 'off') < 0.05;
            anova_p(iCell) = anova1(all_resp_plaid_trials(iCell,:), all_plaid_labels, 'off') < 0.05;
        end

        % --- gDSI (same as exptAnalysis lines 216-217) ---
        fprintf('  Computing gDSI...\n')
        gDSI_grat_struct  = getGDSI(avg_resp(:, 1:nStimDir), stimDirs);
        gDSI_plaid_struct = getGDSI(avg_resp(:, nStimDir+1:end), stimDirs);

        nDS_both = sum(gDSI_grat_struct.gDSI > gDSI_threshold & gDSI_plaid_struct.gDSI > gDSI_threshold);
        fprintf('  gDSI > %.1f: grat %d/%d, plaid %d/%d, both %d/%d\n', gDSI_threshold, ...
            length(gDSI_grat_struct.DS_ind), nCells, ...
            length(gDSI_plaid_struct.DS_ind), nCells, ...
            nDS_both, nCells)

        % --- Peak DSI (for comparison with gDSI) ---
        DSI = zeros(1, nCells);
        for iCell = 1:nCells
            [max_val, max_ind] = max(avg_resp(iCell, 1:nStimDir));
            null_ind = max_ind + (nStimDir/2);
            if null_ind > nStimDir; null_ind = null_ind - nStimDir; end
            min_val = avg_resp(iCell, null_ind);
            if min_val < 0; min_val = 0; end
            DSI(iCell) = (max_val - min_val) / (max_val + min_val);
        end

        % --- Angular distance to VS/IOC predictions ---
        delta_sess = gDSI_grat_struct.prefDir_deg - gDSI_plaid_struct.prefDir_deg;
        delta_sess = mod(delta_sess + 180, 360) - 180;
        vs_angle_sess  = -shift_VS_deg;
        ioc_angle_sess = -shift_IOC_deg;
        dist_VS_sess  = min(abs(delta_sess - vs_angle_sess), 360 - abs(delta_sess - vs_angle_sess));
        dist_IOC_sess = min(abs(delta_sess - ioc_angle_sess), 360 - abs(delta_sess - ioc_angle_sess));
        dist_C2_sess  = min(abs(delta_sess - c2_angle_sess), 360 - abs(delta_sess - c2_angle_sess));

        % --- VS/IOC-resolved Zp/Zc (per-session) ---
        % Reshape flat avg_resp into avg_resp_dir format for getZpZcStruct_VSIOC
        % avg_resp is nCells x 2*nStimDir; need nCells x nStimDir x 1 x 2 x 1
        sess_avg_resp_dir = zeros(nCells, nStimDir, 1, 2, 1);
        sess_avg_resp_dir(:,:,1,1,1) = avg_resp(:, 1:nStimDir);
        sess_avg_resp_dir(:,:,1,2,1) = avg_resp(:, nStimDir+1:end);

        fprintf('  Computing Zp/Zc VSIOC...\n')
        sess_ZpZcVSIOC = getZpZcStruct_VSIOC(sess_avg_resp_dir, alpha_VS, alpha_IOC, 'whole_cell', plaid_offset);
        fprintf('  Zp/Zc VSIOC: degenerate=%d, VS-PDS=%d, IOC-PDS=%d, CDS=%d\n', ...
            sess_ZpZcVSIOC.is_degenerate, length(sess_ZpZcVSIOC.ind_VS_PDS), ...
            length(sess_ZpZcVSIOC.ind_IOC_PDS), length(sess_ZpZcVSIOC.ind_CDS))

        % --- Store session metadata ---
        sessions(idx).mouse           = mouse;
        sessions(idx).date            = date;
        sessions(idx).run_str         = run_str;
        sessions(idx).nCells          = nCells;
        sessions(idx).cell_range      = [cell_offset+1, cell_offset+nCells];
        sessions(idx).speed_stim      = speed_stim;
        sessions(idx).speed_mask      = speed_mask;
        sessions(idx).plaid_offset_deg = plaid_offset;
        sessions(idx).alpha_VS        = alpha_VS;
        sessions(idx).alpha_IOC       = alpha_IOC;
        sessions(idx).shift_VS_deg    = shift_VS_deg;
        sessions(idx).shift_IOC_deg   = shift_IOC_deg;
        sessions(idx).c2_angle_deltapref = c2_angle_sess;
        sessions(idx).shift_C2_deg    = shift_C2_deg_sess;
        sessions(idx).ZpZcVSIOC_degenerate = sess_ZpZcVSIOC.is_degenerate;

        if ~sess_ZpZcVSIOC.is_degenerate
            any_nondegen = true;
        end

        % --- Concatenate ---
        all_resp_mean    = [all_resp_mean; avg_resp];
        all_resp_sem     = [all_resp_sem; avg_resp_se];
        all_gDSI_grat    = [all_gDSI_grat, gDSI_grat_struct.gDSI];
        all_gDSI_plaid   = [all_gDSI_plaid, gDSI_plaid_struct.gDSI];
        all_prefDir_grat = [all_prefDir_grat, gDSI_grat_struct.prefDir_deg];
        all_prefDir_plaid = [all_prefDir_plaid, gDSI_plaid_struct.prefDir_deg];
        all_resp_flag    = [all_resp_flag, resp_flag];
        all_anova_grat   = [all_anova_grat, anova_g];
        all_anova_plaid  = [all_anova_plaid, anova_p];
        all_session_id   = [all_session_id, idx*ones(1, nCells)];
        all_Zp_VS        = [all_Zp_VS, sess_ZpZcVSIOC.Zp_VS(1,:)];
        all_Zp_IOC       = [all_Zp_IOC, sess_ZpZcVSIOC.Zp_IOC(1,:)];
        all_Zc_VSIOC     = [all_Zc_VSIOC, sess_ZpZcVSIOC.Zc(1,:)];
        all_DSI          = [all_DSI, DSI];
        all_dist_VS      = [all_dist_VS, dist_VS_sess];
        all_dist_IOC     = [all_dist_IOC, dist_IOC_sess];
        all_dist_C2      = [all_dist_C2, dist_C2_sess];
        all_resp_any     = [all_resp_any, resp_any_flag];
        cell_offset = cell_offset + nCells;
    end

    %% Build summary struct & save
    fprintf('\n=== Building summary struct (plaid angle = %d deg) ===\n', this_angle)
    fprintf('Total cells: %d across %d sessions\n', cell_offset, length(angle_expts))

    S.stim_dirs        = ref_stimDirs;
    S.plaid_offset_deg = sessions(1).plaid_offset_deg;
    S.resp_mean        = all_resp_mean;
    S.resp_sem         = all_resp_sem;
    S.gDSI_grat        = all_gDSI_grat;
    S.gDSI_plaid       = all_gDSI_plaid;
    S.prefDir_grat_deg = all_prefDir_grat;
    S.prefDir_plaid_deg = all_prefDir_plaid;
    S.resp_to_grat     = logical(all_resp_flag);
    S.anova_grat       = logical(all_anova_grat);
    S.anova_plaid      = logical(all_anova_plaid);
    S.DS_grat          = all_gDSI_grat > gDSI_threshold;
    S.DS_plaid         = all_gDSI_plaid > gDSI_threshold;
    S.DS_both          = S.DS_grat & S.DS_plaid;
    S.session_id       = all_session_id;
    S.sessions         = sessions;
    S.Zp_VS            = all_Zp_VS;
    S.Zp_IOC           = all_Zp_IOC;
    S.Zc_VSIOC         = all_Zc_VSIOC;
    S.any_nondegen_VSIOC = any_nondegen;
    S.DSI_grat    = all_DSI;
    S.dist_VS     = all_dist_VS;
    S.dist_IOC    = all_dist_IOC;
    S.dist_C2     = all_dist_C2;
    S.resp_to_any = logical(all_resp_any);

    outFile = fullfile(outDir, 'aggregate_crossOriSummary_gDSI.mat');
    save(outFile, '-struct', 'S');
    fprintf('Saved: %s\n', outFile)

    %% Group sessions by speed condition (for VS/IOC prediction lines)
    fprintf('\n=== Plotting population figures (plaid angle = %d deg) ===\n', this_angle)
    nSessions = length(S.sessions);
    nTotal    = size(S.resp_mean, 1);

    speeds = zeros(nSessions, 2);
    for i = 1:nSessions
        speeds(i,:) = [S.sessions(i).speed_stim, S.sessions(i).speed_mask];
    end
    [unique_speeds, cond_first_idx, ~] = unique(speeds, 'rows');
    nConditions = size(unique_speeds, 1);

    %% Figure 1: gDSI Scatter (grating vs plaid pref dir)
    % Based on exptAnalysis_DC_gDSI.m lines 823-845

    figure('Visible', 'off', 'Position', [0 0 1080 1080]);
    hold on

    % VS/IOC/C2 prediction lines per speed condition
    h_vs_lines  = gobjects(nConditions, 1);
    h_ioc_lines = gobjects(nConditions, 1);
    vs_colors  = [0 0.7 0; 0 0.4 0];   % green shades
    ioc_colors = [0.8 0 0.8; 0.5 0 0.5]; % magenta shades
    for k = 1:nConditions
        si = cond_first_idx(k);
        sVS  = S.sessions(si).shift_VS_deg;
        sIOC = S.sessions(si).shift_IOC_deg;
        vc = vs_colors(min(k, size(vs_colors,1)), :);
        ic = ioc_colors(min(k, size(ioc_colors,1)), :);
        h_vs_lines(k)  = plot([0 360], [0+sVS 360+sVS],   '--', 'Color', vc, 'LineWidth', 1);
        plot([0 360], [0+sVS-360 360+sVS-360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sVS+360 360+sVS+360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        h_ioc_lines(k) = plot([0 360], [0+sIOC 360+sIOC], '--', 'Color', ic, 'LineWidth', 1);
        plot([0 360], [0+sIOC-360 360+sIOC-360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sIOC+360 360+sIOC+360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
    end

    % Scatter ALL cells, color by session
    session_colors = lines(nSessions);
    h_scatter = gobjects(nSessions, 1);
    for i = 1:nSessions
        mask = S.session_id == i;
        h_scatter(i) = scatter(S.prefDir_grat_deg(mask), S.prefDir_plaid_deg(mask), ...
            30, session_colors(i,:), 'LineWidth', 0.8);
    end

    xlabel('Grating gDSI pref dir (deg)')
    ylabel('Plaid gDSI pref dir (deg)')
    xlim([0 360]); ylim([0 360]); axis square

    % Build legend
    leg_handles = gobjects(0,1);
    leg_labels  = {};
    for k = 1:nConditions
        leg_handles = [leg_handles; h_vs_lines(k); h_ioc_lines(k)];
        leg_labels  = [leg_labels, ...
            {sprintf('VS (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}, ...
            {sprintf('IOC (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}];
    end
    for i = 1:nSessions
        leg_handles = [leg_handles; h_scatter(i)];
        leg_labels  = [leg_labels, {[S.sessions(i).mouse ' ' S.sessions(i).date]}];
    end
    legend(leg_handles, leg_labels, 'Location', 'best', 'FontSize', 7)

    sgtitle({sprintf('Grating vs Plaid preferred direction (gDSI, %d\\circ plaid)', this_angle), ...
        ['n=' num2str(nTotal) ' cells, ' num2str(nSessions) ' sessions']}, 'FontSize', 10)

    print(fullfile(outDir, 'aggregate_GratVsPlaidPrefDir_gDSI_summary.png'), '-dpng', '-r300')
    fprintf('Saved scatter: aggregate_GratVsPlaidPrefDir_gDSI_summary.png\n')

    %% Figure 1b: gDSI Scatter — DS cells only (gDSI > threshold both)
    figure('Visible', 'off', 'Position', [0 0 1080 1080]);
    hold on

    % VS/IOC lines (same as Figure 1)
    h_vs_lines  = gobjects(nConditions, 1);
    h_ioc_lines = gobjects(nConditions, 1);
    for k = 1:nConditions
        si = cond_first_idx(k);
        sVS  = S.sessions(si).shift_VS_deg;
        sIOC = S.sessions(si).shift_IOC_deg;
        vc = vs_colors(min(k, size(vs_colors,1)), :);
        ic = ioc_colors(min(k, size(ioc_colors,1)), :);
        h_vs_lines(k)  = plot([0 360], [0+sVS 360+sVS],   '--', 'Color', vc, 'LineWidth', 1);
        plot([0 360], [0+sVS-360 360+sVS-360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sVS+360 360+sVS+360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        h_ioc_lines(k) = plot([0 360], [0+sIOC 360+sIOC], '--', 'Color', ic, 'LineWidth', 1);
        plot([0 360], [0+sIOC-360 360+sIOC-360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sIOC+360 360+sIOC+360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
    end

    % Scatter DS cells only, color by session
    h_scatter = gobjects(nSessions, 1);
    for i = 1:nSessions
        mask = S.session_id == i & S.DS_both;
        h_scatter(i) = scatter(S.prefDir_grat_deg(mask), S.prefDir_plaid_deg(mask), ...
            30, session_colors(i,:), 'LineWidth', 0.8);
    end

    xlabel('Grating gDSI pref dir (deg)')
    ylabel('Plaid gDSI pref dir (deg)')
    xlim([0 360]); ylim([0 360]); axis square

    % Legend (same structure as Figure 1)
    leg_handles = gobjects(0,1);
    leg_labels  = {};
    for k = 1:nConditions
        leg_handles = [leg_handles; h_vs_lines(k); h_ioc_lines(k)];
        leg_labels  = [leg_labels, ...
            {sprintf('VS (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}, ...
            {sprintf('IOC (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}];
    end
    for i = 1:nSessions
        leg_handles = [leg_handles; h_scatter(i)];
        leg_labels  = [leg_labels, {[S.sessions(i).mouse ' ' S.sessions(i).date]}];
    end
    legend(leg_handles, leg_labels, 'Location', 'best', 'FontSize', 7)

    sgtitle({sprintf('Grating vs Plaid pref dir (gDSI > %.1f both, %d\\circ plaid)', gDSI_threshold, this_angle), ...
        ['n=' num2str(sum(S.DS_both)) '/' num2str(nTotal) ' cells']}, 'FontSize', 10)

    print(fullfile(outDir, 'aggregate_GratVsPlaidPrefDir_gDSI_DS_summary.png'), '-dpng', '-r300')
    fprintf('Saved scatter: aggregate_GratVsPlaidPrefDir_gDSI_DS_summary.png\n')

    %% Figure 1c: gDSI Scatter — colored by gDSI (grating)
    figure('Visible', 'off', 'Position', [0 0 1080 1080]);
    hold on

    h_vs_lines  = gobjects(nConditions, 1);
    h_ioc_lines = gobjects(nConditions, 1);
    for k = 1:nConditions
        si = cond_first_idx(k);
        sVS  = S.sessions(si).shift_VS_deg;
        sIOC = S.sessions(si).shift_IOC_deg;
        vc = vs_colors(min(k, size(vs_colors,1)), :);
        ic = ioc_colors(min(k, size(ioc_colors,1)), :);
        h_vs_lines(k)  = plot([0 360], [0+sVS 360+sVS],   '--', 'Color', vc, 'LineWidth', 1);
        plot([0 360], [0+sVS-360 360+sVS-360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sVS+360 360+sVS+360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        h_ioc_lines(k) = plot([0 360], [0+sIOC 360+sIOC], '--', 'Color', ic, 'LineWidth', 1);
        plot([0 360], [0+sIOC-360 360+sIOC-360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sIOC+360 360+sIOC+360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
    end

    scatter(S.prefDir_grat_deg, S.prefDir_plaid_deg, 30, S.gDSI_grat, 'filled')
    colormap(parula)
    cb = colorbar; cb.Label.String = 'gDSI (grating)';
    set(gca, 'CLim', [0 1])

    xlabel('Grating gDSI pref dir (deg)')
    ylabel('Plaid gDSI pref dir (deg)')
    xlim([0 360]); ylim([0 360]); axis square

    leg_handles = gobjects(0,1);
    leg_labels  = {};
    for k = 1:nConditions
        leg_handles = [leg_handles; h_vs_lines(k); h_ioc_lines(k)];
        leg_labels  = [leg_labels, ...
            {sprintf('VS (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}, ...
            {sprintf('IOC (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}];
    end
    legend(leg_handles, leg_labels, 'Location', 'best', 'FontSize', 7)

    sgtitle({sprintf('Grating vs Plaid pref dir — colored by gDSI (grating), %d\\circ plaid', this_angle), ...
        ['n=' num2str(nTotal) ' cells, ' num2str(nSessions) ' sessions']}, 'FontSize', 10)

    print(fullfile(outDir, 'aggregate_GratVsPlaidPrefDir_gDSI_colorGrat.png'), '-dpng', '-r300')
    fprintf('Saved scatter: aggregate_GratVsPlaidPrefDir_gDSI_colorGrat.png\n')

    %% Figure 1d: gDSI Scatter — colored by gDSI (plaid)
    figure('Visible', 'off', 'Position', [0 0 1080 1080]);
    hold on

    h_vs_lines  = gobjects(nConditions, 1);
    h_ioc_lines = gobjects(nConditions, 1);
    for k = 1:nConditions
        si = cond_first_idx(k);
        sVS  = S.sessions(si).shift_VS_deg;
        sIOC = S.sessions(si).shift_IOC_deg;
        vc = vs_colors(min(k, size(vs_colors,1)), :);
        ic = ioc_colors(min(k, size(ioc_colors,1)), :);
        h_vs_lines(k)  = plot([0 360], [0+sVS 360+sVS],   '--', 'Color', vc, 'LineWidth', 1);
        plot([0 360], [0+sVS-360 360+sVS-360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sVS+360 360+sVS+360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        h_ioc_lines(k) = plot([0 360], [0+sIOC 360+sIOC], '--', 'Color', ic, 'LineWidth', 1);
        plot([0 360], [0+sIOC-360 360+sIOC-360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sIOC+360 360+sIOC+360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
    end

    scatter(S.prefDir_grat_deg, S.prefDir_plaid_deg, 30, S.gDSI_plaid, 'filled')
    colormap(parula)
    cb = colorbar; cb.Label.String = 'gDSI (plaid)';
    set(gca, 'CLim', [0 1])

    xlabel('Grating gDSI pref dir (deg)')
    ylabel('Plaid gDSI pref dir (deg)')
    xlim([0 360]); ylim([0 360]); axis square

    leg_handles = gobjects(0,1);
    leg_labels  = {};
    for k = 1:nConditions
        leg_handles = [leg_handles; h_vs_lines(k); h_ioc_lines(k)];
        leg_labels  = [leg_labels, ...
            {sprintf('VS (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}, ...
            {sprintf('IOC (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}];
    end
    legend(leg_handles, leg_labels, 'Location', 'best', 'FontSize', 7)

    sgtitle({sprintf('Grating vs Plaid pref dir — colored by gDSI (plaid), %d\\circ plaid', this_angle), ...
        ['n=' num2str(nTotal) ' cells, ' num2str(nSessions) ' sessions']}, 'FontSize', 10)

    print(fullfile(outDir, 'aggregate_GratVsPlaidPrefDir_gDSI_colorPlaid.png'), '-dpng', '-r300')
    fprintf('Saved scatter: aggregate_GratVsPlaidPrefDir_gDSI_colorPlaid.png\n')

    %% Figure 1e: gDSI Scatter — DS cells, colored by gDSI (grating)
    figure('Visible', 'off', 'Position', [0 0 1080 1080]);
    hold on

    h_vs_lines  = gobjects(nConditions, 1);
    h_ioc_lines = gobjects(nConditions, 1);
    for k = 1:nConditions
        si = cond_first_idx(k);
        sVS  = S.sessions(si).shift_VS_deg;
        sIOC = S.sessions(si).shift_IOC_deg;
        vc = vs_colors(min(k, size(vs_colors,1)), :);
        ic = ioc_colors(min(k, size(ioc_colors,1)), :);
        h_vs_lines(k)  = plot([0 360], [0+sVS 360+sVS],   '--', 'Color', vc, 'LineWidth', 1);
        plot([0 360], [0+sVS-360 360+sVS-360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sVS+360 360+sVS+360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        h_ioc_lines(k) = plot([0 360], [0+sIOC 360+sIOC], '--', 'Color', ic, 'LineWidth', 1);
        plot([0 360], [0+sIOC-360 360+sIOC-360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sIOC+360 360+sIOC+360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
    end

    scatter(S.prefDir_grat_deg(S.DS_both), S.prefDir_plaid_deg(S.DS_both), 30, S.gDSI_grat(S.DS_both), 'filled')
    colormap(parula)
    cb = colorbar; cb.Label.String = 'gDSI (grating)';
    set(gca, 'CLim', [0 1])

    xlabel('Grating gDSI pref dir (deg)')
    ylabel('Plaid gDSI pref dir (deg)')
    xlim([0 360]); ylim([0 360]); axis square

    leg_handles = gobjects(0,1);
    leg_labels  = {};
    for k = 1:nConditions
        leg_handles = [leg_handles; h_vs_lines(k); h_ioc_lines(k)];
        leg_labels  = [leg_labels, ...
            {sprintf('VS (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}, ...
            {sprintf('IOC (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}];
    end
    legend(leg_handles, leg_labels, 'Location', 'best', 'FontSize', 7)

    sgtitle({sprintf('Grating vs Plaid pref dir — colored by gDSI (grating), %d\\circ plaid', this_angle), ...
        ['gDSI > ' sprintf('%.1f', gDSI_threshold) ' both, n=' num2str(sum(S.DS_both)) '/' num2str(nTotal)]}, 'FontSize', 10)

    print(fullfile(outDir, 'aggregate_GratVsPlaidPrefDir_gDSI_DS_colorGrat.png'), '-dpng', '-r300')
    fprintf('Saved scatter: aggregate_GratVsPlaidPrefDir_gDSI_DS_colorGrat.png\n')

    %% Figure 1f: gDSI Scatter — DS cells, colored by gDSI (plaid)
    figure('Visible', 'off', 'Position', [0 0 1080 1080]);
    hold on

    h_vs_lines  = gobjects(nConditions, 1);
    h_ioc_lines = gobjects(nConditions, 1);
    for k = 1:nConditions
        si = cond_first_idx(k);
        sVS  = S.sessions(si).shift_VS_deg;
        sIOC = S.sessions(si).shift_IOC_deg;
        vc = vs_colors(min(k, size(vs_colors,1)), :);
        ic = ioc_colors(min(k, size(ioc_colors,1)), :);
        h_vs_lines(k)  = plot([0 360], [0+sVS 360+sVS],   '--', 'Color', vc, 'LineWidth', 1);
        plot([0 360], [0+sVS-360 360+sVS-360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sVS+360 360+sVS+360], '--', 'Color', vc, 'LineWidth', 1, 'HandleVisibility', 'off');
        h_ioc_lines(k) = plot([0 360], [0+sIOC 360+sIOC], '--', 'Color', ic, 'LineWidth', 1);
        plot([0 360], [0+sIOC-360 360+sIOC-360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
        plot([0 360], [0+sIOC+360 360+sIOC+360], '--', 'Color', ic, 'LineWidth', 1, 'HandleVisibility', 'off');
    end

    scatter(S.prefDir_grat_deg(S.DS_both), S.prefDir_plaid_deg(S.DS_both), 30, S.gDSI_plaid(S.DS_both), 'filled')
    colormap(parula)
    cb = colorbar; cb.Label.String = 'gDSI (plaid)';
    set(gca, 'CLim', [0 1])

    xlabel('Grating gDSI pref dir (deg)')
    ylabel('Plaid gDSI pref dir (deg)')
    xlim([0 360]); ylim([0 360]); axis square

    leg_handles = gobjects(0,1);
    leg_labels  = {};
    for k = 1:nConditions
        leg_handles = [leg_handles; h_vs_lines(k); h_ioc_lines(k)];
        leg_labels  = [leg_labels, ...
            {sprintf('VS (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}, ...
            {sprintf('IOC (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}];
    end
    legend(leg_handles, leg_labels, 'Location', 'best', 'FontSize', 7)

    sgtitle({sprintf('Grating vs Plaid pref dir — colored by gDSI (plaid), %d\\circ plaid', this_angle), ...
        ['gDSI > ' sprintf('%.1f', gDSI_threshold) ' both, n=' num2str(sum(S.DS_both)) '/' num2str(nTotal)]}, 'FontSize', 10)

    print(fullfile(outDir, 'aggregate_GratVsPlaidPrefDir_gDSI_DS_colorPlaid.png'), '-dpng', '-r300')
    fprintf('Saved scatter: aggregate_GratVsPlaidPrefDir_gDSI_DS_colorPlaid.png\n')

    %% Figure 2: gDSI Polar Histogram (delta pref dir) with 4-category classification
    % Based on exptAnalysis_DC_gDSI.m, updated with C2 prediction

    half_window_deg = (360 / nStimDir) / 4;  % quarter stimulus step (±5.625° for 16 dirs)

    % Classify per-session (each session has its own VS/IOC/C2 angles)
    DS_both_idx = find(S.DS_both);
    delta = S.prefDir_grat_deg(DS_both_idx) - S.prefDir_plaid_deg(DS_both_idx);
    delta = mod(delta + 180, 360) - 180;  % wrap to [-180, 180]

    % Per-session 4-category classification, then combine
    all_clf_labels = repmat({''}, 1, length(delta));
    agg_n_VS = 0; agg_n_IOC = 0; agg_n_C1 = 0; agg_n_C2 = 0; agg_n_patt = 0; agg_n_unc = 0;
    agg_degenerate_pairs = {};
    for i = 1:nSessions
        sess_mask_DS = S.session_id(DS_both_idx) == i;
        if ~any(sess_mask_DS), continue; end

        sess_vs_angle  = -S.sessions(i).shift_VS_deg;
        sess_ioc_angle = -S.sessions(i).shift_IOC_deg;
        sess_c2_angle  = S.sessions(i).c2_angle_deltapref;
        sess_delta = delta(sess_mask_DS);

        sess_clf = classifyDeltaPrefDir(sess_delta, sess_vs_angle, sess_ioc_angle, sess_c2_angle, half_window_deg);
        sess_indices = find(sess_mask_DS);
        for j = 1:length(sess_indices)
            all_clf_labels{sess_indices(j)} = sess_clf.labels{j};
        end
        agg_n_VS  = agg_n_VS  + sess_clf.n_VS;
        agg_n_IOC = agg_n_IOC + sess_clf.n_IOC;
        agg_n_C1  = agg_n_C1  + sess_clf.n_C1;
        agg_n_C2  = agg_n_C2  + sess_clf.n_C2;
        agg_n_patt = agg_n_patt + sess_clf.n_pattern;
        agg_n_unc = agg_n_unc + sess_clf.n_unclass;
        if ~isempty(sess_clf.degenerate_pairs) && isempty(agg_degenerate_pairs)
            agg_degenerate_pairs = sess_clf.degenerate_pairs;
        end
    end

    figure('Visible', 'off', 'Position', [0 0 720 720]);
    polarhistogram(deg2rad(delta), deg2rad(-180:11.25:180), ...
        'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 0.6, 'EdgeColor', 'k')
    hold on

    r_max = rlim;
    r_max = r_max(2);

    % VS/IOC/C2 radial lines per speed condition
    leg_h = [];  leg_l = {};
    for k = 1:nConditions
        si = cond_first_idx(k);
        sVS  = S.sessions(si).shift_VS_deg;
        sIOC = S.sessions(si).shift_IOC_deg;
        sC2_angle = S.sessions(si).c2_angle_deltapref;
        vs_angle  = -sVS;
        ioc_angle = -sIOC;

        h1 = polarplot(deg2rad([vs_angle vs_angle]),   [0 r_max], 'b--', 'LineWidth', 1.5);
        h2 = polarplot(deg2rad([ioc_angle ioc_angle]), [0 r_max], 'r--', 'LineWidth', 1.5);
        if nConditions > 1
            leg_h = [leg_h, h1, h2];
            leg_l = [leg_l, {sprintf('VS (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}, ...
                            {sprintf('IOC (v_s=%.0f, v_m=%.0f)', unique_speeds(k,1), unique_speeds(k,2))}];
        else
            leg_h = [leg_h, h1, h2];
            leg_l = [leg_l, {'VS'}, {'IOC'}];
        end
    end
    legend(leg_h, leg_l, 'Location', 'best', 'FontSize', 8)

    % Colored scatter ring at r_max * 1.08
    r_ring = r_max * 1.08;
    % Find cells containing each label (labels may be merged like 'VS/C1')
    idx_VS   = cellfun(@(l) contains(l, 'VS') && ~contains(l, 'IOC'), all_clf_labels);
    idx_IOC  = cellfun(@(l) contains(l, 'IOC') && ~contains(l, 'VS'), all_clf_labels);
    idx_C1   = cellfun(@(l) contains(l, 'C1') && ~contains(l, 'VS') && ~contains(l, 'IOC'), all_clf_labels);
    idx_patt = cellfun(@(l) contains(l, 'VS') && contains(l, 'IOC'), all_clf_labels);
    idx_unc  = strcmp(all_clf_labels, 'unclassified');

    if any(idx_VS)
        polarplot(deg2rad(delta(idx_VS)), r_ring*ones(1,sum(idx_VS)), ...
            'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.2 0.4 0.8], 'MarkerEdgeColor', 'none')
    end
    if any(idx_IOC)
        polarplot(deg2rad(delta(idx_IOC)), r_ring*ones(1,sum(idx_IOC)), ...
            'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.8 0.2 0.2], 'MarkerEdgeColor', 'none')
    end
    if any(idx_C1)
        polarplot(deg2rad(delta(idx_C1)), r_ring*ones(1,sum(idx_C1)), ...
            'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.4 0.4 0.4], 'MarkerEdgeColor', 'none')
    end
    if any(idx_patt)
        polarplot(deg2rad(delta(idx_patt)), r_ring*ones(1,sum(idx_patt)), ...
            'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.5 0.3 0.7], 'MarkerEdgeColor', 'none')
    end
    if any(idx_unc)
        polarplot(deg2rad(delta(idx_unc)), r_ring*ones(1,sum(idx_unc)), ...
            'o', 'MarkerSize', 5, 'MarkerFaceColor', [0.75 0.75 0.75], 'MarkerEdgeColor', 'none')
    end

    % Build classification string for title
    clf_str = sprintf('VS:%d  IOC:%d  C1:%d', agg_n_VS, agg_n_IOC, agg_n_C1);
    if agg_n_patt > 0
        clf_str = [clf_str sprintf('  Pattern(merged):%d', agg_n_patt)];
    end
    clf_str = [clf_str sprintf('  Unclass:%d', agg_n_unc)];
    if ~isempty(agg_degenerate_pairs)
        degen_strs = cellfun(@(p) [p{1} '/' p{2}], agg_degenerate_pairs, 'UniformOutput', false);
        clf_str = [clf_str '  [near-degen: ' strjoin(degen_strs, ', ') ']'];
    end

    sgtitle({sprintf('\\DeltaPrefDir (grating - plaid), gDSI > %.1f both, %d\\circ plaid', gDSI_threshold, this_angle), ...
        ['n=' num2str(sum(S.DS_both)) '/' num2str(nTotal) ', window=\pm' num2str(half_window_deg) '\circ'], ...
        clf_str}, 'FontSize', 9)

    fprintf('  Delta pref dir classification: %s\n', clf_str)

    % Save classification to struct
    S.clf_labels = all_clf_labels;
    S.clf_DS_idx = DS_both_idx;

    print(fullfile(outDir, 'aggregate_DeltaPrefDirPolar_gDSI_summary.png'), '-dpng', '-r300')
    fprintf('Saved polar: aggregate_DeltaPrefDirPolar_gDSI_summary.png\n')

    % Store for cross-angle KDE overlay (plotted after loop)
    kde_delta_by_angle{ai} = delta;
    kde_angle_labels(ai)   = this_angle;
    kde_n_DS(ai)           = length(delta);

    %% Figure 2b: Aggregate Cell Classification Pie Chart & Bar
    fprintf('=== Plotting aggregate cell classification pie chart (plaid angle = %d deg) ===\n', this_angle)

    pie_labels = {'VS', 'IOC', 'C1', 'Unclassified'};
    pie_counts = [agg_n_VS, agg_n_IOC, agg_n_C1, agg_n_unc + agg_n_C2];
    pie_colors = [0.2 0.4 0.8; 0.8 0.2 0.2; 0.4 0.4 0.4; 0.75 0.75 0.75];
    if agg_n_patt > 0
        pie_labels = [pie_labels, {'VS/IOC (merged)'}];
        pie_counts = [pie_counts, agg_n_patt];
        pie_colors = [pie_colors; 0.5 0.3 0.7];
    end

    nz = pie_counts > 0;
    pie_labels_nz = pie_labels(nz);
    pie_counts_nz = pie_counts(nz);
    pie_colors_nz = pie_colors(nz, :);
    n_total_clf = sum(pie_counts_nz);

    figure('Visible', 'off', 'Position', [0 0 1200 500]);
    subplot(1, 2, 1)
    if ~isempty(pie_counts_nz) && any(pie_counts_nz > 0)
        h_pie = pie(pie_counts_nz);
        for pp = 1:length(pie_counts_nz)
            h_pie(2*pp-1).FaceColor = pie_colors_nz(pp, :);
        end
        legend(pie_labels_nz, 'Location', 'bestoutside', 'FontSize', 8)
    end
    title(sprintf('Classification (n=%d, window=\\pm%d\\circ)', n_total_clf, half_window_deg))

    subplot(1, 2, 2)
    bh = barh(pie_counts_nz, 'FaceColor', 'flat');
    bh.CData = pie_colors_nz;
    yticks(1:length(pie_counts_nz))
    yticklabels(pie_labels_nz)
    xlabel('Number of cells')
    for pp = 1:length(pie_counts_nz)
        text(pie_counts_nz(pp) + max(pie_counts_nz)*0.03, pp, ...
            sprintf('%d (%.1f%%)', pie_counts_nz(pp), 100*pie_counts_nz(pp)/n_total_clf), ...
            'VerticalAlignment', 'middle', 'FontSize', 9)
    end
    xlim([0 max(pie_counts_nz)*1.3])
    title('Cell counts by category')

    sgtitle({sprintf('Aggregate Cell Classification — %d\\circ plaid, %d sessions', this_angle, nSessions), ...
        sprintf('n=%d DS cells, half-window=%d\\circ', n_total_clf, half_window_deg)}, 'FontSize', 10)

    print(fullfile(outDir, 'aggregate_CellClassification_gDSI.png'), '-dpng', '-r300')
    fprintf('Saved: aggregate_CellClassification_gDSI.png\n')

    %% Aggregate Circular V-test
    fprintf('=== Aggregate circular V-test (plaid angle = %d deg) ===\n', this_angle)

    alpha_rad_agg = deg2rad(delta);
    n_vtest_agg = length(alpha_rad_agg);

    % Use first session's predictions (all sessions in same angle group share same offset)
    ref_vs_angle  = -S.sessions(1).shift_VS_deg;
    ref_ioc_angle = -S.sessions(1).shift_IOC_deg;
    ref_c2_angle  = S.sessions(1).c2_angle_deltapref;

    vtest_agg = struct();
    vtest_names = {'VS', 'IOC', 'C1', 'C2'};
    vtest_angles = [ref_vs_angle, ref_ioc_angle, 0, ref_c2_angle];

    for vv = 1:4
        mu_rad = deg2rad(vtest_angles(vv));
        C_bar = mean(cos(alpha_rad_agg - mu_rad));
        V_stat = n_vtest_agg * C_bar;
        u_stat = V_stat * sqrt(2 / n_vtest_agg);
        p_val = exp(sqrt(1 + 4*n_vtest_agg + 4*(n_vtest_agg^2 - V_stat^2)) - (1 + 2*n_vtest_agg));
        p_val = min(p_val, 1);

        vtest_agg.(vtest_names{vv}).V = V_stat;
        vtest_agg.(vtest_names{vv}).u = u_stat;
        vtest_agg.(vtest_names{vv}).p = p_val;
        vtest_agg.(vtest_names{vv}).mu_deg = vtest_angles(vv);

        fprintf('  V-test toward %s (%+.1f deg): V=%.1f, u=%.2f, p=%.2e\n', ...
            vtest_names{vv}, vtest_angles(vv), V_stat, u_stat, p_val)
    end

    V_vals_agg = [vtest_agg.VS.V, vtest_agg.IOC.V, vtest_agg.C1.V, vtest_agg.C2.V];
    [~, winner_idx_agg] = max(V_vals_agg);
    vtest_winner_agg = vtest_names{winner_idx_agg};
    fprintf('  V-test winner: %s (V=%.1f, p=%.2e)\n', ...
        vtest_winner_agg, V_vals_agg(winner_idx_agg), vtest_agg.(vtest_winner_agg).p)

    % Save V-test results to struct
    S.vtest_VS  = vtest_agg.VS;
    S.vtest_IOC = vtest_agg.IOC;
    S.vtest_C1  = vtest_agg.C1;
    S.vtest_C2  = vtest_agg.C2;

    %% Aggregate gDSI Threshold Sweep
    fprintf('=== Aggregate gDSI threshold sweep (plaid angle = %d deg) ===\n', this_angle)

    sweep_thresholds = [0.2, 0.3, 0.4, 0.5];
    nThresh = length(sweep_thresholds);
    sweep_n_cells  = zeros(1, nThresh);
    sweep_n_VS     = zeros(1, nThresh);
    sweep_n_IOC    = zeros(1, nThresh);
    sweep_n_C1     = zeros(1, nThresh);
    sweep_n_C2     = zeros(1, nThresh);
    sweep_n_unc    = zeros(1, nThresh);
    sweep_V_VS     = zeros(1, nThresh);
    sweep_V_IOC    = zeros(1, nThresh);

    for tt = 1:nThresh
        thr = sweep_thresholds(tt);
        ind_sweep = find(S.gDSI_grat > thr & S.gDSI_plaid > thr);
        sweep_n_cells(tt) = length(ind_sweep);

        if sweep_n_cells(tt) < 3, continue; end

        dp_sweep = S.prefDir_grat_deg(ind_sweep) - S.prefDir_plaid_deg(ind_sweep);
        dp_sweep = mod(dp_sweep + 180, 360) - 180;

        % Per-session classification at this threshold
        sw_n_VS = 0; sw_n_IOC = 0; sw_n_C1 = 0; sw_n_C2 = 0; sw_n_unc = 0;
        for i = 1:nSessions
            sess_mask = S.session_id(ind_sweep) == i;
            if ~any(sess_mask), continue; end
            sess_clf = classifyDeltaPrefDir(dp_sweep(sess_mask), ...
                -S.sessions(i).shift_VS_deg, -S.sessions(i).shift_IOC_deg, ...
                S.sessions(i).c2_angle_deltapref, half_window_deg);
            sw_n_VS  = sw_n_VS  + sess_clf.n_VS;
            sw_n_IOC = sw_n_IOC + sess_clf.n_IOC;
            sw_n_C1  = sw_n_C1  + sess_clf.n_C1;
            sw_n_C2  = sw_n_C2  + sess_clf.n_C2;
            sw_n_unc = sw_n_unc + sess_clf.n_unclass + sess_clf.n_pattern;
        end
        sweep_n_VS(tt)  = sw_n_VS;
        sweep_n_IOC(tt) = sw_n_IOC;
        sweep_n_C1(tt)  = sw_n_C1;
        sweep_n_C2(tt)  = sw_n_C2;
        sweep_n_unc(tt) = sw_n_unc;

        % V-test at each threshold
        a_sweep = deg2rad(dp_sweep);
        n_sw = length(a_sweep);
        sweep_V_VS(tt)  = n_sw * mean(cos(a_sweep - deg2rad(ref_vs_angle)));
        sweep_V_IOC(tt) = n_sw * mean(cos(a_sweep - deg2rad(ref_ioc_angle)));

        fprintf('  gDSI > %.1f: n=%d, VS:%d IOC:%d C1:%d C2:%d unc:%d, V_VS=%.1f V_IOC=%.1f\n', ...
            thr, sweep_n_cells(tt), sw_n_VS, sw_n_IOC, sw_n_C1, sw_n_C2, sw_n_unc, ...
            sweep_V_VS(tt), sweep_V_IOC(tt))
    end

    figure('Visible', 'off', 'Position', [0 0 1200 900]);

    % Top-left: Stacked bar of proportions
    subplot(2, 2, 1)
    prop_data = [sweep_n_VS; sweep_n_IOC; sweep_n_C1; sweep_n_unc + sweep_n_C2]';
    bar_h = bar(prop_data, 'stacked');
    bar_h(1).FaceColor = [0.2 0.4 0.8];
    bar_h(2).FaceColor = [0.8 0.2 0.2];
    bar_h(3).FaceColor = [0.4 0.4 0.4];
    bar_h(4).FaceColor = [0.75 0.75 0.75];
    xticks(1:nThresh)
    xticklabels(arrayfun(@(x) sprintf('%.1f', x), sweep_thresholds, 'UniformOutput', false))
    xlabel('gDSI threshold')
    ylabel('Number of cells')
    legend({'VS', 'IOC', 'C1', 'Unclass'}, 'Location', 'bestoutside', 'FontSize', 7)
    title('Classification by threshold')

    % Top-right: V-statistic for VS and IOC
    subplot(2, 2, 2)
    plot(sweep_thresholds, sweep_V_VS, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b')
    hold on
    plot(sweep_thresholds, sweep_V_IOC, 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r')
    xlabel('gDSI threshold')
    ylabel('V-statistic')
    legend({'V toward VS', 'V toward IOC'}, 'Location', 'best', 'FontSize', 8)
    title('V-test strength vs threshold')

    % Bottom-left: n_cells remaining
    subplot(2, 2, 3)
    bar(sweep_thresholds, sweep_n_cells, 0.5, 'FaceColor', [0.3 0.3 0.6])
    xlabel('gDSI threshold')
    ylabel('n cells (gDSI both > thr)')
    title('Cells surviving threshold')
    for tt = 1:nThresh
        text(sweep_thresholds(tt), sweep_n_cells(tt) + max(sweep_n_cells)*0.03, ...
            num2str(sweep_n_cells(tt)), 'HorizontalAlignment', 'center', 'FontSize', 9)
    end

    % Bottom-right: Overlaid histograms at each threshold
    subplot(2, 2, 4)
    hold on
    colors_sweep = [0.2 0.4 0.8; 0.4 0.7 0.4; 0.9 0.6 0.1; 0.8 0.2 0.2];
    for tt = 1:nThresh
        thr = sweep_thresholds(tt);
        ind_sw = find(S.gDSI_grat > thr & S.gDSI_plaid > thr);
        if length(ind_sw) < 3, continue; end
        dp_sw = S.prefDir_grat_deg(ind_sw) - S.prefDir_plaid_deg(ind_sw);
        dp_sw = mod(dp_sw + 180, 360) - 180;
        [counts, edges] = histcounts(dp_sw, -180:11.25:180);
        centers = (edges(1:end-1) + edges(2:end)) / 2;
        counts_norm = counts / sum(counts);
        plot(centers, counts_norm, '-', 'Color', colors_sweep(tt,:), 'LineWidth', 1.5)
    end
    xline(ref_vs_angle, 'b--', 'LineWidth', 1)
    xline(ref_ioc_angle, 'r--', 'LineWidth', 1)
    xline(0, 'k--', 'LineWidth', 0.5)
    xlabel('\DeltaPrefDir (deg)')
    ylabel('Normalized density')
    legend([arrayfun(@(x) sprintf('gDSI>%.1f', x), sweep_thresholds, 'UniformOutput', false), ...
        {'VS', 'IOC', 'C1'}], 'Location', 'best', 'FontSize', 6)
    title('Distribution shape vs threshold')
    xlim([-180 180])

    sgtitle({sprintf('Aggregate gDSI Threshold Sweep — %d\\circ plaid, %d sessions', this_angle, nSessions), ...
        sprintf('half-window=%d\\circ', half_window_deg)}, 'FontSize', 10)

    print(fullfile(outDir, 'aggregate_gDSI_ThresholdSweep.png'), '-dpng', '-r300')
    fprintf('Saved: aggregate_gDSI_ThresholdSweep.png\n')

    % Re-save struct with all new fields
    save(outFile, '-struct', 'S');
    fprintf('Re-saved with classification + V-test: %s\n', outFile)

    %% Figure 3: Zp/Zc VSIOC scatter (only when non-degenerate sessions exist)
    if any_nondegen
        fprintf('\n=== Plotting Zp/Zc VSIOC figures (plaid angle = %d deg) ===\n', this_angle)

        % Classify aggregate Zp/Zc
        threshold = 1.28;
        agg_VS_PDS  = (S.Zp_VS > threshold) & (S.Zp_VS - S.Zc_VSIOC > threshold);
        agg_IOC_PDS = (S.Zp_IOC > threshold) & (S.Zp_IOC - S.Zc_VSIOC > threshold);
        agg_CDS     = (S.Zc_VSIOC > threshold) & (S.Zc_VSIOC - max(S.Zp_VS, S.Zp_IOC) > threshold);

        % Resolve ties
        both_pat = agg_VS_PDS & agg_IOC_PDS;
        vs_wins  = both_pat & (S.Zp_VS >= S.Zp_IOC);
        ioc_wins = both_pat & (S.Zp_IOC > S.Zp_VS);
        agg_VS_PDS  = (agg_VS_PDS & ~both_pat) | vs_wins;
        agg_IOC_PDS = (agg_IOC_PDS & ~both_pat) | ioc_wins;

        n_VS_PDS  = sum(agg_VS_PDS);
        n_IOC_PDS = sum(agg_IOC_PDS);
        n_CDS_agg = sum(agg_CDS);

        % Figure 3a: Zp_VS vs Zp_IOC scatter
        figure('Visible', 'off', 'Position', [0 0 720 720]);
        hold on
        plot([-4 8], [-4 8], 'k--', 'LineWidth', 0.5)
        plot([-4 8], [1.28 1.28], 'k:', 'LineWidth', 0.5)
        plot([1.28 1.28], [-4 8], 'k:', 'LineWidth', 0.5)

        scatter(S.Zp_VS, S.Zp_IOC, 20, [0.75 0.75 0.75], 'filled', 'MarkerFaceAlpha', 0.5)
        if n_VS_PDS > 0
            scatter(S.Zp_VS(agg_VS_PDS), S.Zp_IOC(agg_VS_PDS), 30, [0.2 0.4 0.8], 'filled')
        end
        if n_IOC_PDS > 0
            scatter(S.Zp_VS(agg_IOC_PDS), S.Zp_IOC(agg_IOC_PDS), 30, [0.8 0.2 0.2], 'filled')
        end
        if n_CDS_agg > 0
            scatter(S.Zp_VS(agg_CDS), S.Zp_IOC(agg_CDS), 30, [0.4 0.4 0.4], 'filled')
        end

        xlabel('Zp_{VS}'); ylabel('Zp_{IOC}')
        xlim([-4 8]); ylim([-4 8]); axis square
        legend({'y=x', '', '', 'all cells', ...
            sprintf('VS-pattern (n=%d)', n_VS_PDS), ...
            sprintf('IOC-pattern (n=%d)', n_IOC_PDS), ...
            sprintf('Component (n=%d)', n_CDS_agg)}, ...
            'Location', 'best', 'FontSize', 7)
        sgtitle({sprintf('Zp_{VS} vs Zp_{IOC} (aggregate, %d\\circ plaid)', this_angle), ...
            sprintf('VS-PDS:%d, IOC-PDS:%d, CDS:%d / %d cells', n_VS_PDS, n_IOC_PDS, n_CDS_agg, nTotal)}, 'FontSize', 10)

        print(fullfile(outDir, 'aggregate_ZpVS_vs_ZpIOC_gDSI.png'), '-dpng', '-r300')
        fprintf('Saved scatter: aggregate_ZpVS_vs_ZpIOC_gDSI.png\n')

        % Figure 3b: Side-by-side Zp_VS vs Zc and Zp_IOC vs Zc
        figure('Visible', 'off', 'Position', [0 0 1440 720]);
        subplot(1,2,1)
            hold on
            scatter(S.Zc_VSIOC, S.Zp_VS, 20, [0.75 0.75 0.75], 'filled', 'MarkerFaceAlpha', 0.5)
            if n_VS_PDS > 0
                scatter(S.Zc_VSIOC(agg_VS_PDS), S.Zp_VS(agg_VS_PDS), 30, [0.2 0.4 0.8], 'filled')
            end
            if n_CDS_agg > 0
                scatter(S.Zc_VSIOC(agg_CDS), S.Zp_VS(agg_CDS), 30, [0.4 0.4 0.4], 'filled')
            end
            plotZcZpBorders
            xlabel('Zc'); ylabel('Zp_{VS}')
            xlim([-4 8]); ylim([-4 8]); axis square
            title(sprintf('VS-pattern: %d PDS, %d CDS', n_VS_PDS, n_CDS_agg))
        subplot(1,2,2)
            hold on
            scatter(S.Zc_VSIOC, S.Zp_IOC, 20, [0.75 0.75 0.75], 'filled', 'MarkerFaceAlpha', 0.5)
            if n_IOC_PDS > 0
                scatter(S.Zc_VSIOC(agg_IOC_PDS), S.Zp_IOC(agg_IOC_PDS), 30, [0.8 0.2 0.2], 'filled')
            end
            if n_CDS_agg > 0
                scatter(S.Zc_VSIOC(agg_CDS), S.Zp_IOC(agg_CDS), 30, [0.4 0.4 0.4], 'filled')
            end
            plotZcZpBorders
            xlabel('Zc'); ylabel('Zp_{IOC}')
            xlim([-4 8]); ylim([-4 8]); axis square
            title(sprintf('IOC-pattern: %d PDS, %d CDS', n_IOC_PDS, n_CDS_agg))

        sgtitle(sprintf('Zp vs Zc (VS/IOC resolved, aggregate, %d\\circ plaid)', this_angle), 'FontSize', 10)
        print(fullfile(outDir, 'aggregate_ZpZc_VSIOC_gDSI.png'), '-dpng', '-r300')
        fprintf('Saved scatter: aggregate_ZpZc_VSIOC_gDSI.png\n')
    else
        fprintf('  Skipping Zp/Zc VSIOC figures (all sessions degenerate for this angle group)\n')
    end

    %% Figure 4: Cell Filtering Summary Funnel
    fprintf('\n=== Plotting cell filtering funnel (plaid angle = %d deg) ===\n', this_angle)

    funnel_labels = {'Total imaged', 'Resp (any)', 'Resp (grating)', ...
        'ANOVA dir-sel', 'Resp + ANOVA', 'Resp + ANOVA + gDSI', 'gDSI both'};
    funnel_counts = [nTotal, sum(S.resp_to_any), sum(S.resp_to_grat), ...
        sum(S.anova_grat), sum(S.resp_to_grat & S.anova_grat), ...
        sum(S.resp_to_grat & S.anova_grat & S.DS_grat), sum(S.DS_both)];

    figure('Visible', 'off', 'Position', [0 0 900 500]);
    bh = barh(flip(funnel_counts), 'FaceColor', 'flat');
    colors = [linspace(0.85, 0.1, 7); linspace(0.85, 0.3, 7); linspace(0.85, 0.7, 7)]';
    bh.CData = flip(colors);
    yticks(1:7)
    yticklabels(flip(funnel_labels))
    xlabel('Number of cells')
    for i = 1:7
        text(funnel_counts(8-i) + max(funnel_counts)*0.02, i, ...
            sprintf('%d (%.0f%%)', funnel_counts(8-i), 100*funnel_counts(8-i)/nTotal), ...
            'VerticalAlignment', 'middle', 'FontSize', 9)
    end
    title(sprintf('Cell Filtering Summary — %d\\circ plaid, %d sessions', this_angle, nSessions), 'FontSize', 11)
    xlim([0 max(funnel_counts)*1.25])

    print(fullfile(outDir, 'aggregate_CellFilteringSummary_gDSI.png'), '-dpng', '-r300')
    fprintf('Saved funnel: aggregate_CellFilteringSummary_gDSI.png\n')

    %% Figure 5: gDSI Histograms
    fprintf('=== Plotting gDSI histograms (plaid angle = %d deg) ===\n', this_angle)

    figure('Visible', 'off', 'Position', [0 0 1080 720]);
    subplot(2,1,1)
        histogram(S.gDSI_grat)
        hold on
        xline(gDSI_threshold, 'r--', 'LineWidth', 1.5)
        xlabel('gDSI')
        ylabel('# cells')
        title(sprintf('Grating gDSI — %d/%d cells > %.1f', sum(S.DS_grat), nTotal, gDSI_threshold))
    subplot(2,1,2)
        histogram(S.gDSI_plaid)
        hold on
        xline(gDSI_threshold, 'r--', 'LineWidth', 1.5)
        xlabel('gDSI')
        ylabel('# cells')
        title(sprintf('Plaid gDSI — %d/%d cells > %.1f', sum(S.DS_plaid), nTotal, gDSI_threshold))
    sgtitle(sprintf('gDSI distributions — %d\\circ plaid, %d sessions, n=%d', ...
        this_angle, nSessions, nTotal), 'FontSize', 10)

    print(fullfile(outDir, 'aggregate_gDSI_histograms_gDSI.png'), '-dpng', '-r300')
    fprintf('Saved histograms: aggregate_gDSI_histograms_gDSI.png\n')

    %% Figure 6: DSI vs gDSI Scatter with VS/IOC Distance Gradient
    fprintf('=== Plotting DSI vs gDSI gradient scatter (plaid angle = %d deg) ===\n', this_angle)

    ind_scatter = S.resp_to_grat & S.anova_grat & isfinite(S.DSI_grat);
    n_scatter = sum(ind_scatter);

    figure('Visible', 'off', 'Position', [0 0 1440 720]);

    % Left: colored by VS distance
    subplot(1,2,1)
    scatter(S.gDSI_grat(ind_scatter), S.DSI_grat(ind_scatter), 30, S.dist_VS(ind_scatter), 'filled')
    hold on
    plot([0 1], [0 1], 'k--', 'LineWidth', 0.5)
    colormap(subplot(1,2,1), flipud(jet))
    cb = colorbar; cb.Label.String = 'Angular distance to VS (deg)';
    clim([0 180])
    xlim([0 1]); ylim([0 1]); axis square
    xlabel('gDSI (grating)'); ylabel('Peak DSI (grating)')
    title(sprintf('VS distance  (n=%d)', n_scatter))

    % Right: colored by IOC distance
    subplot(1,2,2)
    scatter(S.gDSI_grat(ind_scatter), S.DSI_grat(ind_scatter), 30, S.dist_IOC(ind_scatter), 'filled')
    hold on
    plot([0 1], [0 1], 'k--', 'LineWidth', 0.5)
    colormap(subplot(1,2,2), flipud(jet))
    cb = colorbar; cb.Label.String = 'Angular distance to IOC (deg)';
    clim([0 180])
    xlim([0 1]); ylim([0 1]); axis square
    xlabel('gDSI (grating)'); ylabel('Peak DSI (grating)')
    title(sprintf('IOC distance  (n=%d)', n_scatter))

    sgtitle({sprintf('DSI vs gDSI (grating) — %d\\circ plaid, %d sessions', this_angle, nSessions), ...
        sprintf('resp + ANOVA, n=%d/%d', n_scatter, nTotal)}, 'FontSize', 10)

    print(fullfile(outDir, 'aggregate_DSI_vs_gDSI_gradient_gDSI.png'), '-dpng', '-r300')
    fprintf('Saved scatter: aggregate_DSI_vs_gDSI_gradient_gDSI.png\n')

end  % angle loop

%% Cross-angle KDE overlay: von Mises kernel density of delta pref dir
% Von Mises kernel is the circular analog of Gaussian KDE — correctly
% handles ±180° wraparound, no bin-size parameter, single bandwidth κ.
% κ = 10 gives ~18° smoothing width, matching 22.5° stimulus spacing.
fprintf('\n=== Plotting cross-angle von Mises KDE overlay ===\n')

kdeOutDir = fullfile(base, 'Analysis\2P\aggregate');
if ~exist(kdeOutDir, 'dir'), mkdir(kdeOutDir); end

kappa = 10;                              % concentration parameter (tunable)
theta_eval = linspace(-180, 180, 361);   % 1° resolution
theta_eval_rad = deg2rad(theta_eval);
I0_kappa = besseli(0, kappa);

valid_mask = kde_n_DS >= 3;
nValid = sum(valid_mask);

if nValid == 0
    fprintf('  No angle groups with >= 3 DS cells — skipping KDE overlay.\n')
else
    kde_colors = lines(length(unique_angles));

    figure('Visible', 'off', 'Position', [0 0 900 500]);
    hold on

    for ai = 1:length(unique_angles)
        if ~valid_mask(ai), continue; end

        d_rad = deg2rad(kde_delta_by_angle{ai});
        n = length(d_rad);

        % Vectorized von Mises KDE
        diff_mat = theta_eval_rad(:) - d_rad(:)';     % 361 × n
        kde_vals = sum(exp(kappa * cos(diff_mat)), 2)' / (n * 2 * pi * I0_kappa);
        kde_vals_deg = kde_vals * (pi / 180);          % density per degree

        plot(theta_eval, kde_vals_deg, '-', ...
            'Color', kde_colors(ai,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('%d\\circ plaid (n=%d)', kde_angle_labels(ai), n))

        fprintf('  %d deg: n=%d DS cells, peak=%.4f /deg\n', ...
            kde_angle_labels(ai), n, max(kde_vals_deg))
    end

    xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    xlabel('\DeltaPrefDir (grating - plaid, deg)')
    ylabel('Density (per degree)')
    xlim([-180 180])
    legend('Location', 'best', 'FontSize', 9)
    title(sprintf('Von Mises KDE of \\DeltaPrefDir (\\kappa=%d, gDSI > %.1f)', ...
        kappa, gDSI_threshold), 'FontSize', 10)

    print(fullfile(kdeOutDir, 'aggregate_DeltaPrefDir_KDE_overlay.png'), '-dpng', '-r300')
    fprintf('Saved: aggregate_DeltaPrefDir_KDE_overlay.png\n')
end

%% Cross-angle smoothed histogram overlay of delta pref dir
% Simpler alternative to the von Mises KDE: bin into equal-width bins,
% normalize to fraction of cells, then Gaussian-smooth with circular
% padding to handle ±180° wraparound.
fprintf('\n=== Plotting cross-angle smoothed histogram overlay ===\n')

hist_bin_width = 1;                                     % degrees per bin
hist_bin_edges = -180:hist_bin_width:180;               % 360 bins
hist_bin_centers = hist_bin_edges(1:end-1) + hist_bin_width/2;
hist_smooth_window = 30;                                % bins (~30 deg Gaussian window)
hist_colors = lines(length(unique_angles) + 1);         % offset by 1 for distinct palette
hist_colors = hist_colors(2:end, :);                    % skip first color (used by KDE)

if nValid == 0
    fprintf('  No angle groups with >= 3 DS cells — skipping smoothed histogram.\n')
else
    figure('Visible', 'off', 'Position', [0 0 900 500]);
    hold on

    for ai = 1:length(unique_angles)
        if ~valid_mask(ai), continue; end

        delta_deg = kde_delta_by_angle{ai};
        n = length(delta_deg);

        % Bin into histogram and normalize to fraction
        counts = histcounts(delta_deg, hist_bin_edges);
        frac = counts / n;

        % Circular padding for wraparound-safe smoothing
        pad = hist_smooth_window;
        frac_padded = [frac(end-pad+1:end), frac, frac(1:pad)];
        frac_smooth = smoothdata(frac_padded, 'gaussian', hist_smooth_window);
        frac_smooth = frac_smooth(pad+1:end-pad);       % trim back

        plot(hist_bin_centers, frac_smooth, '-', ...
            'Color', hist_colors(ai,:), 'LineWidth', 2, ...
            'DisplayName', sprintf('%d\\circ plaid (n=%d)', kde_angle_labels(ai), n))

        fprintf('  %d deg: n=%d DS cells, peak fraction=%.3f\n', ...
            kde_angle_labels(ai), n, max(frac_smooth))
    end

    xline(0, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    xlabel('\DeltaPrefDir (grating - plaid, deg)')
    ylabel('Fraction of cells')
    xlim([-180 180])
    legend('Location', 'best', 'FontSize', 9)
    title(sprintf('Smoothed \\DeltaPrefDir histogram (%d\\circ bins, gDSI > %.1f)', ...
        hist_bin_width, gDSI_threshold), 'FontSize', 10)

    print(fullfile(kdeOutDir, 'aggregate_DeltaPrefDir_SmoothedHist_overlay.png'), '-dpng', '-r300')
    fprintf('Saved: aggregate_DeltaPrefDir_SmoothedHist_overlay.png\n')
end

fprintf('\nDone.\n')
