function [ZpZcVSIOC] = getZpZcStruct_VSIOC(avg_resp_dir, alpha_VS, alpha_IOC, data_type, plaid_offset_deg)
% getZpZcStruct_VSIOC  VS/IOC-resolved partial correlations
%
%   ZpZcVSIOC = getZpZcStruct_VSIOC(avg_resp_dir, alpha_VS, alpha_IOC)
%   ZpZcVSIOC = getZpZcStruct_VSIOC(avg_resp_dir, alpha_VS, alpha_IOC, data_type, plaid_offset_deg)
%
% Extends the standard Zp/Zc framework to distinguish VS-pattern from
% IOC-pattern cells using paired 2-predictor partial correlations.
%
% Input:
%   avg_resp_dir    - nCells x nDir x nMaskPhase x 2(grat/plaid) x 2(mean/std)
%   alpha_VS        - scalar (radians), VS angular shift
%   alpha_IOC       - scalar (radians), IOC angular shift
%   data_type       - 'whole_cell' (default) or 'alignedTestDir'
%   plaid_offset_deg - plaid offset in degrees (default 45)
%
% Output struct fields:
%   Zp_VS, Zp_IOC, Zc  - nMaskPhase x nCells partial correlation Z-scores
%   Rp_VS, Rp_IOC, Rc_A, Rc_B - raw partial correlations
%   VS_pattern_pred, IOC_pattern_pred, component_pred - prediction arrays
%   ind_VS_PDS, ind_IOC_PDS, ind_CDS - classification indices (across phases)
%   ind_VS_PDS_byphase, ind_IOC_PDS_byphase, ind_CDS_byphase - per phase
%   is_degenerate   - true if alpha_VS ≈ alpha_IOC (fell back to standard)
%   standard_ZpZc   - standard getZpZcStruct output (when degenerate)

    if nargin < 4 || isempty(data_type)
        data_type = 'whole_cell';
    end
    if nargin < 5 || isempty(plaid_offset_deg)
        plaid_offset_deg = 45;
    end

    nCells    = size(avg_resp_dir, 1);
    nStimDir  = size(avg_resp_dir, 2);
    nMaskPhas = size(avg_resp_dir, 3);
    int = 360 / nStimDir;  % degrees per bin

    % --- Degeneracy check ---
    sep_deg = abs(rad2deg(alpha_VS) - rad2deg(alpha_IOC));
    sep_deg = min(sep_deg, 360 - sep_deg);
    is_degenerate = sep_deg < 1.0;

    if is_degenerate
        fprintf('  getZpZcStruct_VSIOC: VS ≈ IOC (sep=%.2f deg) — falling back to standard getZpZcStruct\n', sep_deg)
        standard = getZpZcStruct(avg_resp_dir, data_type, plaid_offset_deg);

        ZpZcVSIOC.Zp_VS  = standard.Zp;
        ZpZcVSIOC.Zp_IOC = standard.Zp;  % identical when degenerate
        ZpZcVSIOC.Zc     = standard.Zc;
        ZpZcVSIOC.Rp_VS  = standard.Rp;
        ZpZcVSIOC.Rp_IOC = standard.Rp;
        ZpZcVSIOC.Rc_A   = standard.Rc;
        ZpZcVSIOC.Rc_B   = standard.Rc;
        ZpZcVSIOC.VS_pattern_pred  = standard.pattern_pred;
        ZpZcVSIOC.IOC_pattern_pred = standard.pattern_pred;
        ZpZcVSIOC.component_pred   = standard.component_pred;
        ZpZcVSIOC.ind_VS_PDS       = standard.PDSind_all;
        ZpZcVSIOC.ind_IOC_PDS      = standard.PDSind_all;
        ZpZcVSIOC.ind_CDS          = standard.CDSind_all;
        ZpZcVSIOC.ind_VS_PDS_byphase  = standard.PDSind_byphase;
        ZpZcVSIOC.ind_IOC_PDS_byphase = standard.PDSind_byphase;
        ZpZcVSIOC.ind_CDS_byphase     = standard.CDSind_byphase;
        ZpZcVSIOC.is_degenerate = true;
        ZpZcVSIOC.standard_ZpZc = standard;
        return
    end

    % --- Build prediction arrays ---
    grating_resp = avg_resp_dir(:,:,1,1,1);  % nCells x nDir

    % VS pattern prediction: shift grating by alpha_VS
    shift_VS_bins = rad2deg(alpha_VS) / int;
    VS_pattern_pred = circshift_interp(grating_resp, -shift_VS_bins);

    % IOC pattern prediction: shift grating by alpha_IOC
    shift_IOC_bins = rad2deg(alpha_IOC) / int;
    IOC_pattern_pred = circshift_interp(grating_resp, -shift_IOC_bins);

    % Component prediction (same as standard getZpZcStruct)
    half_offset = plaid_offset_deg / 2;
    if strcmp(data_type, 'alignedTestDir')
        component = grating_resp + circshift(grating_resp, -plaid_offset_deg/int, 2);
    elseif strcmp(data_type, 'whole_cell')
        component = circshift(grating_resp, +half_offset/int, 2) + ...
                    circshift(grating_resp, -half_offset/int, 2);
    end

    % --- Analysis A: VS-pattern vs component ---
    comp_corr_A      = zeros(nMaskPhas, nCells);
    patt_corr_VS     = zeros(nMaskPhas, nCells);
    comp_patt_corr_A = zeros(nMaskPhas, nCells);

    % --- Analysis B: IOC-pattern vs component ---
    comp_corr_B       = zeros(nMaskPhas, nCells);
    patt_corr_IOC     = zeros(nMaskPhas, nCells);
    comp_patt_corr_B  = zeros(nMaskPhas, nCells);

    for iCell = 1:nCells
        for ip = 1:nMaskPhas
            plaid_resp = avg_resp_dir(iCell,:,ip,2,1);

            % Analysis A (VS vs component)
            comp_corr_A(ip,iCell)      = triu2vec(corrcoef(plaid_resp, component(iCell,:)));
            patt_corr_VS(ip,iCell)     = triu2vec(corrcoef(plaid_resp, VS_pattern_pred(iCell,:)));
            comp_patt_corr_A(ip,iCell) = triu2vec(corrcoef(component(iCell,:), VS_pattern_pred(iCell,:)));

            % Analysis B (IOC vs component)
            comp_corr_B(ip,iCell)       = triu2vec(corrcoef(plaid_resp, component(iCell,:)));
            patt_corr_IOC(ip,iCell)     = triu2vec(corrcoef(plaid_resp, IOC_pattern_pred(iCell,:)));
            comp_patt_corr_B(ip,iCell)  = triu2vec(corrcoef(component(iCell,:), IOC_pattern_pred(iCell,:)));
        end
    end

    % --- Partial correlations (Analysis A: VS) ---
    Rp_VS = (patt_corr_VS - comp_corr_A .* comp_patt_corr_A) ./ ...
            sqrt((1 - comp_corr_A.^2) .* (1 - comp_patt_corr_A.^2));
    Rc_A  = (comp_corr_A - patt_corr_VS .* comp_patt_corr_A) ./ ...
            sqrt((1 - patt_corr_VS.^2) .* (1 - comp_patt_corr_A.^2));

    % --- Partial correlations (Analysis B: IOC) ---
    Rp_IOC = (patt_corr_IOC - comp_corr_B .* comp_patt_corr_B) ./ ...
             sqrt((1 - comp_corr_B.^2) .* (1 - comp_patt_corr_B.^2));
    Rc_B   = (comp_corr_B - patt_corr_IOC .* comp_patt_corr_B) ./ ...
             sqrt((1 - patt_corr_IOC.^2) .* (1 - comp_patt_corr_B.^2));

    % --- Fisher Z transformation ---
    denom = sqrt(1 / (nStimDir - 3));
    Zp_VS  = (0.5 .* log((1 + Rp_VS)  ./ (1 - Rp_VS)))  ./ denom;
    Zp_IOC = (0.5 .* log((1 + Rp_IOC) ./ (1 - Rp_IOC))) ./ denom;
    Zc_A   = (0.5 .* log((1 + Rc_A)   ./ (1 - Rc_A)))   ./ denom;
    Zc_B   = (0.5 .* log((1 + Rc_B)   ./ (1 - Rc_B)))   ./ denom;

    % Combined Zc = mean of both analyses
    Zc = (Zc_A + Zc_B) / 2;

    % --- Classification ---
    threshold = 1.28;
    ind_VS_PDS_byphase  = cell(1, nMaskPhas);
    ind_IOC_PDS_byphase = cell(1, nMaskPhas);
    ind_CDS_byphase     = cell(1, nMaskPhas);

    for ip = 1:nMaskPhas
        % VS-pattern: Zp_VS > threshold AND Zp_VS - Zc > threshold
        vs_cells  = (Zp_VS(ip,:) > threshold) & (Zp_VS(ip,:) - Zc(ip,:) > threshold);
        % IOC-pattern: Zp_IOC > threshold AND Zp_IOC - Zc > threshold
        ioc_cells = (Zp_IOC(ip,:) > threshold) & (Zp_IOC(ip,:) - Zc(ip,:) > threshold);
        % Component: Zc > threshold AND Zc - max(Zp_VS, Zp_IOC) > threshold
        cds_cells = (Zc(ip,:) > threshold) & (Zc(ip,:) - max(Zp_VS(ip,:), Zp_IOC(ip,:)) > threshold);

        % Resolve ties: if both VS and IOC satisfied, winner = higher Zp
        both_pat = vs_cells & ioc_cells;
        vs_wins  = both_pat & (Zp_VS(ip,:) >= Zp_IOC(ip,:));
        ioc_wins = both_pat & (Zp_IOC(ip,:) > Zp_VS(ip,:));
        vs_cells  = (vs_cells & ~both_pat) | vs_wins;
        ioc_cells = (ioc_cells & ~both_pat) | ioc_wins;

        ind_VS_PDS_byphase{ip}  = find(vs_cells);
        ind_IOC_PDS_byphase{ip} = find(ioc_cells);
        ind_CDS_byphase{ip}     = find(cds_cells);
    end

    ind_VS_PDS  = unique([ind_VS_PDS_byphase{:}]);
    ind_IOC_PDS = unique([ind_IOC_PDS_byphase{:}]);
    ind_CDS     = unique([ind_CDS_byphase{:}]);

    % --- Pack output ---
    ZpZcVSIOC.Zp_VS    = Zp_VS;
    ZpZcVSIOC.Zp_IOC   = Zp_IOC;
    ZpZcVSIOC.Zc        = Zc;
    ZpZcVSIOC.Rp_VS    = Rp_VS;
    ZpZcVSIOC.Rp_IOC   = Rp_IOC;
    ZpZcVSIOC.Rc_A     = Rc_A;
    ZpZcVSIOC.Rc_B     = Rc_B;
    ZpZcVSIOC.VS_pattern_pred  = VS_pattern_pred;
    ZpZcVSIOC.IOC_pattern_pred = IOC_pattern_pred;
    ZpZcVSIOC.component_pred   = component;
    ZpZcVSIOC.ind_VS_PDS       = ind_VS_PDS;
    ZpZcVSIOC.ind_IOC_PDS      = ind_IOC_PDS;
    ZpZcVSIOC.ind_CDS          = ind_CDS;
    ZpZcVSIOC.ind_VS_PDS_byphase  = ind_VS_PDS_byphase;
    ZpZcVSIOC.ind_IOC_PDS_byphase = ind_IOC_PDS_byphase;
    ZpZcVSIOC.ind_CDS_byphase     = ind_CDS_byphase;
    ZpZcVSIOC.is_degenerate = false;
    ZpZcVSIOC.standard_ZpZc = [];
end
