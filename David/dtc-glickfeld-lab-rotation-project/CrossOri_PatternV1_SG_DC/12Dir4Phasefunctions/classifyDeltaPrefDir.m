function [clf] = classifyDeltaPrefDir(delta_pref_deg, vs_angle_deg, ioc_angle_deg, c2_angle_deg, half_window_deg)
% classifyDeltaPrefDir  Classify cells on delta pref dir polar histogram
%
%   clf = classifyDeltaPrefDir(delta_pref_deg, vs_angle_deg, ioc_angle_deg)
%   clf = classifyDeltaPrefDir(delta_pref_deg, vs_angle_deg, ioc_angle_deg, half_window_deg)
%       3-category mode (VS / IOC / C1) — backward-compatible
%
%   clf = classifyDeltaPrefDir(delta_pref_deg, vs_angle_deg, ioc_angle_deg, c2_angle_deg, half_window_deg)
%       4-category mode (VS / IOC / C1 / C2) — c2_angle_deg is C2 prediction
%
% Inputs:
%   delta_pref_deg   - 1 x N, wrapped to [-180, 180]
%   vs_angle_deg     - scalar, VS prediction angle (degrees)
%   ioc_angle_deg    - scalar, IOC prediction angle (degrees)
%   c2_angle_deg     - scalar, C2 (mask-tracking) prediction angle (degrees)
%                      Pass [] or omit for 3-category mode.
%   half_window_deg  - scalar, classification window half-width (default 15)
%
% Returns struct clf with fields:
%   ind_VS, ind_IOC, ind_C1, ind_C2 - indices of classified cells
%   ind_comp                         - alias for ind_C1 (backward compat)
%   ind_pattern                      - merged VS+IOC when degenerate (empty otherwise)
%   ind_unclass                      - unclassified cells
%   n_VS, n_IOC, n_C1, n_C2, n_comp, n_pattern, n_unclass - counts
%   dist_VS, dist_IOC, dist_C1, dist_C2 - 1xN angular distance to each prediction
%   is_degenerate                    - true when VS ~ IOC (merged to "pattern")
%   degenerate_pairs                 - cell array of near-degenerate prediction pairs
%   labels                           - 1 x N cell array of labels per cell

    % --- Parse arguments for backward compatibility ---
    % Old signature: (delta, vs, ioc)               → 3-cat, default window
    % Old signature: (delta, vs, ioc, half_window)   → 3-cat, custom window
    % New signature: (delta, vs, ioc, c2, half_window) → 4-cat
    use_c2 = false;
    if nargin < 4
        c2_angle_deg = [];
        half_window_deg = 15;
    elseif nargin == 4
        % 4th arg: could be c2_angle or half_window (old signature)
        % Old code used 5.625 as default — treat scalar < 45 as half_window
        % for backward compat. Angles are typically >= 0 and can be large.
        % Safest: if 4th arg is empty, it's c2=[] with default window.
        if isempty(c2_angle_deg)
            half_window_deg = 15;
        else
            % Old callers pass half_window as 4th arg. Treat as half_window.
            half_window_deg = c2_angle_deg;
            c2_angle_deg = [];
        end
    else  % nargin >= 5
        use_c2 = ~isempty(c2_angle_deg);
        if isempty(half_window_deg)
            half_window_deg = 15;
        end
    end

    N = length(delta_pref_deg);

    % Circular angular distance (degrees, range [0, 180])
    circ_dist_fn = @(a, b) min(abs(a - b), 360 - abs(a - b));

    % --- Compute per-cell distances to all predictions ---
    dist_VS   = circ_dist_fn(delta_pref_deg, vs_angle_deg);
    dist_IOC  = circ_dist_fn(delta_pref_deg, ioc_angle_deg);
    dist_C1   = circ_dist_fn(delta_pref_deg, 0);  % C1 = component tracking test grating = 0 shift
    if use_c2
        dist_C2 = circ_dist_fn(delta_pref_deg, c2_angle_deg);
    else
        dist_C2 = nan(1, N);
    end

    % --- Detect pairwise degeneracy ---
    if use_c2
        pred_names = {'VS', 'IOC', 'C1', 'C2'};
        pred_angles = [vs_angle_deg, ioc_angle_deg, 0, c2_angle_deg];
    else
        pred_names = {'VS', 'IOC', 'C1'};
        pred_angles = [vs_angle_deg, ioc_angle_deg, 0];
    end
    nPred = length(pred_names);

    degenerate_pairs = {};
    degen_matrix = false(nPred);
    for ii = 1:nPred
        for jj = ii+1:nPred
            sep = circ_dist_fn(pred_angles(ii), pred_angles(jj));
            if sep < 2 * half_window_deg
                degenerate_pairs{end+1} = {pred_names{ii}, pred_names{jj}}; %#ok<AGROW>
                degen_matrix(ii, jj) = true;
                degen_matrix(jj, ii) = true;
            end
        end
    end

    % VS-IOC degeneracy (for backward compat "pattern" label)
    vs_ioc_sep = circ_dist_fn(vs_angle_deg, ioc_angle_deg);
    is_degenerate = vs_ioc_sep < 2 * half_window_deg;

    % --- Build merged groups from degeneracy ---
    % Find connected components among degenerate predictions
    visited = false(1, nPred);
    merged_groups = {};  % each entry: indices of predictions that merge
    for ii = 1:nPred
        if visited(ii), continue; end
        group = ii;
        queue = ii;
        visited(ii) = true;
        while ~isempty(queue)
            curr = queue(1); queue(1) = [];
            neighbors = find(degen_matrix(curr, :) & ~visited);
            visited(neighbors) = true;
            group = [group, neighbors]; %#ok<AGROW>
            queue = [queue, neighbors]; %#ok<AGROW>
        end
        merged_groups{end+1} = sort(group); %#ok<AGROW>
    end

    % Build distance matrix: nPred x N
    all_dists = [dist_VS; dist_IOC; dist_C1];
    if use_c2
        all_dists = [all_dists; dist_C2];
    end

    % --- Classify each cell ---
    ind_VS      = [];
    ind_IOC     = [];
    ind_C1      = [];
    ind_C2      = [];
    ind_pattern = [];
    labels      = repmat({''}, 1, N);

    for i = 1:N
        % Find which merged group has the minimum distance
        best_group = 0;
        best_dist  = Inf;
        for g = 1:length(merged_groups)
            group_min = min(all_dists(merged_groups{g}, i));
            if group_min < best_dist
                best_dist = group_min;
                best_group = g;
            end
        end

        if best_dist > half_window_deg
            labels{i} = 'unclassified';
            continue
        end

        group = merged_groups{best_group};
        group_names = pred_names(group);

        if length(group) == 1
            % Single prediction — straightforward
            lbl = group_names{1};
        else
            % Merged group — build combined label
            lbl = strjoin(sort(group_names), '/');
        end

        labels{i} = lbl;

        % Assign to canonical index arrays
        if any(strcmp(group_names, 'VS')) && any(strcmp(group_names, 'IOC'))
            % VS+IOC merged → "pattern" (backward compat)
            ind_pattern = [ind_pattern, i]; %#ok<AGROW>
        end
        if any(strcmp(group_names, 'VS')) && ~any(strcmp(group_names, 'IOC'))
            % VS only (possibly merged with C1 or C2 but not IOC)
            ind_VS = [ind_VS, i]; %#ok<AGROW>
        end
        if any(strcmp(group_names, 'IOC')) && ~any(strcmp(group_names, 'VS'))
            ind_IOC = [ind_IOC, i]; %#ok<AGROW>
        end
        if any(strcmp(group_names, 'C1'))
            ind_C1 = [ind_C1, i]; %#ok<AGROW>
        end
        if use_c2 && any(strcmp(group_names, 'C2'))
            ind_C2 = [ind_C2, i]; %#ok<AGROW>
        end
    end

    ind_unclass = find(strcmp(labels, 'unclassified'));

    % --- Build output struct ---
    clf.ind_VS       = ind_VS;
    clf.ind_IOC      = ind_IOC;
    clf.ind_C1       = ind_C1;
    clf.ind_comp     = ind_C1;        % backward-compat alias
    clf.ind_C2       = ind_C2;
    clf.ind_pattern  = ind_pattern;
    clf.ind_unclass  = ind_unclass;
    clf.n_VS         = length(ind_VS);
    clf.n_IOC        = length(ind_IOC);
    clf.n_C1         = length(ind_C1);
    clf.n_comp       = length(ind_C1); % backward-compat alias
    clf.n_C2         = length(ind_C2);
    clf.n_pattern    = length(ind_pattern);
    clf.n_unclass    = length(ind_unclass);
    clf.dist_VS      = dist_VS;
    clf.dist_IOC     = dist_IOC;
    clf.dist_C1      = dist_C1;
    clf.dist_C2      = dist_C2;
    clf.is_degenerate     = is_degenerate;
    clf.degenerate_pairs  = degenerate_pairs;
    clf.labels       = labels;
end
