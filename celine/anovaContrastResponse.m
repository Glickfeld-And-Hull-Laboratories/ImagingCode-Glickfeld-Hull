function [anova_results, stats_table] = anovaContrastResponse(data1, data2, cell_indices1, cell_indices2, contrasts, sizes, varargin)
% ANOVACONTRASTRESPONSE Perform repeated measures ANOVA on contrast response data
%
% Usage:
%   [anova_results, stats_table] = anovaContrastResponse(data1, data2, cell_indices1, cell_indices2, contrasts, sizes)
%   [anova_results, stats_table] = anovaContrastResponse(data1, data2, cell_indices1, cell_indices2, contrasts, sizes, 'Name', Value, ...)
%
% Inputs:
%   data1 - First dataset (e.g., pref_responses_stat_concat) - pre-drug data
%   data2 - Second dataset - post-drug data (or same as data1 if combined)
%   cell_indices1 - Indices of cells to use for first cell population
%   cell_indices2 - Indices of cells to use for second cell population  
%   contrasts - Array of contrast values
%   sizes - Array of size values
%
% Optional Name-Value Pairs:
%   'Alpha' - Significance level (default: 0.05)
%   'PostHoc' - Whether to perform post-hoc tests (default: true)
%   'Display' - Whether to display results (default: true)
%   'SaveResults' - Path to save results (default: '' - no saving)
%   'CellPopNames' - Names for cell populations {'Pop1', 'Pop2'} (default: {'Population1', 'Population2'})
%
% Outputs:
%   anova_results - Structure containing ANOVA results for each cell population
%   stats_table - Summary table of key statistics

% Parse inputs
p = inputParser;
addRequired(p, 'data1', @iscell);
addRequired(p, 'data2', @iscell);
addRequired(p, 'cell_indices1', @isnumeric);
addRequired(p, 'cell_indices2', @isnumeric);
addRequired(p, 'contrasts', @isnumeric);
addRequired(p, 'sizes', @isnumeric);
addParameter(p, 'Alpha', 0.05, @(x) isnumeric(x) && x > 0 && x < 1);
addParameter(p, 'PostHoc', true, @islogical);
addParameter(p, 'Display', true, @islogical);
addParameter(p, 'SaveResults', '', @ischar);
addParameter(p, 'CellPopNames', {'Population1', 'Population2'}, @(x) iscell(x) && length(x)==2);

parse(p, data1, data2, cell_indices1, cell_indices2, contrasts, sizes, varargin{:});

alpha = p.Results.Alpha;
do_posthoc = p.Results.PostHoc;
display_results = p.Results.Display;
save_path = p.Results.SaveResults;
pop_names = p.Results.CellPopNames;

% Get dimensions
nSizes = length(sizes);
nContrasts = length(contrasts);

% Ask user about size handling if multiple sizes
include_size_effect = false;
if nSizes > 1 && display_results
    fprintf('\nMultiple sizes detected (%d sizes)\n', nSizes);
    response = input('Include Size as a main effect? (y/n): ', 's');
    include_size_effect = strcmpi(response, 'y') || strcmpi(response, 'yes');
    
    if include_size_effect
        fprintf('Will include Size as a main effect in ANOVA\n');
    else
        fprintf('Will run separate ANOVAs for each size\n');
    end
end

% Initialize results
anova_results = struct();
summary_data = {};
row_idx = 1;

if display_results
    fprintf('\n=== Repeated Measures ANOVA Analysis ===\n');
    fprintf('Contrasts: [%s]\n', num2str(contrasts));
    fprintf('Sizes: [%s]\n', num2str(sizes));
    fprintf('Population 1 (%s): %d cells\n', pop_names{1}, length(cell_indices1));
    fprintf('Population 2 (%s): %d cells\n', pop_names{2}, length(cell_indices2));
    fprintf('Alpha level: %.3f\n\n', alpha);
end

% Analyze each cell population
cell_pops = {cell_indices1, cell_indices2};

for pop = 1:2
    current_indices = cell_pops{pop};
    pop_name = pop_names{pop};
    
    if isempty(current_indices)
        if display_results
            fprintf('No cells in %s - skipping\n', pop_name);
        end
        continue;
    end
    
    if include_size_effect && nSizes > 1
        % Single ANOVA with Size as a factor
        field_name = sprintf('%s_AllSizes', strrep(pop_name, ' ', ''));
        
        % Prepare data for repeated measures ANOVA
        rm_data = [];
        cell_id = [];
        
        for iCell = 1:length(current_indices)
            cell_idx = current_indices(iCell);
            cell_data = [];
            
            % Extract data for this cell across all conditions
            for iSize = 1:nSizes
                for iContrast = 1:nContrasts
                    pre_resp = data1{2}(cell_idx, iContrast, iSize);  % Pre-drug (day 2)
                    post_resp = data1{1}(cell_idx, iContrast, iSize); % Post-drug (day 1)
                    
                    if ~isnan(pre_resp) && ~isnan(post_resp)
                        cell_data = [cell_data; pre_resp, post_resp];
                    end
                end
            end
            
            if ~isempty(cell_data)
                rm_data = [rm_data; cell_data];
                cell_id = [cell_id; repmat(iCell, size(cell_data, 1), 1)];
            end
        end
        
        if size(rm_data, 1) < 10
            if display_results
                fprintf('Insufficient data for %s - skipping\n', pop_name);
            end
            continue;
        end
        
        % Create within-subject factors
        contrast_factor = repmat(contrasts', nSizes, 1);
        size_factor = repelem(sizes', nContrasts, 1);
        drug_factor = [1; 2]; % Pre, Post
        
        % Create factor table
        within_design = table(contrast_factor, size_factor, drug_factor(1)*ones(length(contrast_factor),1), ...
            'VariableNames', {'Contrast', 'Size', 'Drug'});
        within_design = [within_design; ...
            table(contrast_factor, size_factor, drug_factor(2)*ones(length(contrast_factor),1), ...
            'VariableNames', {'Contrast', 'Size', 'Drug'})];
        
        % Create variable names for the repeated measures
        var_names = {};
        for iSize = 1:nSizes
            for iContrast = 1:nContrasts
                var_names{end+1} = sprintf('Pre_S%d_C%.2f', iSize, contrasts(iContrast));
                var_names{end+1} = sprintf('Post_S%d_C%.2f', iSize, contrasts(iContrast));
            end
        end
        
        % Reshape data for fitrm
        n_conditions = nSizes * nContrasts * 2;
        reshaped_data = zeros(length(current_indices), n_conditions);
        valid_cells = [];
        
        for iCell = 1:length(current_indices)
            cell_idx = current_indices(iCell);
            col = 1;
            valid_data = true;
            
            for iSize = 1:nSizes
                for iContrast = 1:nContrasts
                    pre_resp = data1{2}(cell_idx, iContrast, iSize);
                    post_resp = data1{1}(cell_idx, iContrast, iSize);
                    
                    if isnan(pre_resp) || isnan(post_resp)
                        valid_data = false;
                        break;
                    end
                    
                    reshaped_data(iCell, col) = pre_resp;
                    reshaped_data(iCell, col+1) = post_resp;
                    col = col + 2;
                end
                if ~valid_data, break; end
            end
            
            if valid_data
                valid_cells = [valid_cells, iCell];
            end
        end
        
        if length(valid_cells) < 5
            if display_results
                fprintf('Insufficient valid cells for %s - skipping\n', pop_name);
            end
            continue;
        end
        
        % Keep only valid cells
        reshaped_data = reshaped_data(valid_cells, :);
        rm_table = array2table(reshaped_data, 'VariableNames', var_names);
        
        try
            % Fit repeated measures model
            rm = fitrm(rm_table, [var_names{1} '-' var_names{end} '~1'], 'WithinDesign', within_design);
            anova_table = ranova(rm, 'WithinModel', 'Contrast*Size*Drug');
            
            % Store results
            anova_results.(field_name).rm_model = rm;
            anova_results.(field_name).anova_table = anova_table;
            anova_results.(field_name).n_cells = length(valid_cells);
            
            % Extract p-values
            contrast_p = anova_table.pValue(strcmp(anova_table.Var1, '(Intercept):Contrast'));
            size_p = anova_table.pValue(strcmp(anova_table.Var1, '(Intercept):Size'));
            drug_p = anova_table.pValue(strcmp(anova_table.Var1, '(Intercept):Drug'));
            contrast_drug_p = anova_table.pValue(strcmp(anova_table.Var1, '(Intercept):Contrast:Drug'));
            size_drug_p = anova_table.pValue(strcmp(anova_table.Var1, '(Intercept):Size:Drug'));
            contrast_size_p = anova_table.pValue(strcmp(anova_table.Var1, '(Intercept):Contrast:Size'));
            three_way_p = anova_table.pValue(strcmp(anova_table.Var1, '(Intercept):Contrast:Size:Drug'));
            
            % Add to summary
            summary_data{row_idx, 1} = pop_name;
            summary_data{row_idx, 2} = 'All';
            summary_data{row_idx, 3} = contrast_p;
            summary_data{row_idx, 4} = drug_p;
            summary_data{row_idx, 5} = contrast_drug_p;
            summary_data{row_idx, 6} = size_p;
            summary_data{row_idx, 7} = length(valid_cells);
            row_idx = row_idx + 1;
            
            if display_results
                fprintf('\n--- %s (All Sizes) ---\n', pop_name);
                fprintf('N cells: %d\n', length(valid_cells));
                fprintf('Contrast effect: p = %.4f %s\n', contrast_p, getSignificanceString(contrast_p, alpha));
                fprintf('Size effect: p = %.4f %s\n', size_p, getSignificanceString(size_p, alpha));
                fprintf('Drug effect: p = %.4f %s\n', drug_p, getSignificanceString(drug_p, alpha));
                fprintf('Contrast × Drug: p = %.4f %s\n', contrast_drug_p, getSignificanceString(contrast_drug_p, alpha));
                fprintf('Size × Drug: p = %.4f %s\n', size_drug_p, getSignificanceString(size_drug_p, alpha));
                fprintf('Contrast × Size: p = %.4f %s\n', contrast_size_p, getSignificanceString(contrast_size_p, alpha));
                fprintf('Three-way interaction: p = %.4f %s\n', three_way_p, getSignificanceString(three_way_p, alpha));
            end
            
        catch ME
            if display_results
                fprintf('ANOVA failed for %s: %s\n', pop_name, ME.message);
            end
        end
        
    else
        % Separate ANOVA for each size
        for iSize = 1:nSizes
            field_name = sprintf('%s_Size%d', strrep(pop_name, ' ', ''), iSize);
            
            % Prepare data for this size
            cell_data = [];
            valid_cells = [];
            
            for iCell = 1:length(current_indices)
                cell_idx = current_indices(iCell);
                pre_data = squeeze(data1{2}(cell_idx, :, iSize));  % Pre-drug
                post_data = squeeze(data1{1}(cell_idx, :, iSize)); % Post-drug
                
                if all(~isnan(pre_data)) && all(~isnan(post_data))
                    % Arrange as [pre_c1, pre_c2, pre_c3, post_c1, post_c2, post_c3]
                    cell_row = [pre_data(:)', post_data(:)'];
                    cell_data = [cell_data; cell_row];
                    valid_cells = [valid_cells, iCell];
                end
            end
            
            if size(cell_data, 1) < 5
                if display_results
                    fprintf('Insufficient data for %s, Size %d - skipping\n', pop_name, sizes(iSize));
                end
                continue;
            end
            
            % Create variable names (must be valid MATLAB variable names)
            var_names = {};
            for iContrast = 1:nContrasts
                var_names{end+1} = sprintf('Pre_C%d', iContrast);
            end
            for iContrast = 1:nContrasts
                var_names{end+1} = sprintf('Post_C%d', iContrast);
            end
            
            % Debug: Check dimensions match
            if size(cell_data, 2) ~= length(var_names)
                if display_results
                    fprintf('Dimension mismatch for %s, Size %d: cell_data has %d columns, var_names has %d elements\n', ...
                        pop_name, sizes(iSize), size(cell_data, 2), length(var_names));
                end
                continue;
            end
            
            % Create within-subject design
            contrast_factor = [contrasts'; contrasts'];
            drug_factor = [ones(nContrasts,1); 2*ones(nContrasts,1)];
            within_design = table(contrast_factor, drug_factor, ...
                'VariableNames', {'Contrast', 'Drug'});
            
            rm_table = array2table(cell_data, 'VariableNames', var_names);
            
            try
                % Fit repeated measures model - use range specification
                formula_str = sprintf('%s-%s~1', var_names{1}, var_names{end});
                rm = fitrm(rm_table, formula_str, 'WithinDesign', within_design);
                anova_table = ranova(rm, 'WithinModel', 'Contrast*Drug');
                
                % Store results
                anova_results.(field_name).rm_model = rm;
                anova_results.(field_name).anova_table = anova_table;
                anova_results.(field_name).n_cells = size(cell_data, 1);
                
                % Extract p-values - check actual column names first
                if display_results
                    fprintf('ANOVA table columns: %s\n', strjoin(anova_table.Properties.VariableNames, ', '));
                end
                
                % Find the correct rows for each effect
                contrast_idx = contains(anova_table.Properties.RowNames, 'Contrast') & ...
                              ~contains(anova_table.Properties.RowNames, 'Drug');
                drug_idx = contains(anova_table.Properties.RowNames, 'Drug') & ...
                          ~contains(anova_table.Properties.RowNames, 'Contrast');
                interaction_idx = contains(anova_table.Properties.RowNames, 'Contrast') & ...
                                 contains(anova_table.Properties.RowNames, 'Drug');
                
                % Extract p-values safely
                if any(contrast_idx)
                    contrast_p = anova_table.pValue(contrast_idx);
                    contrast_p = contrast_p(1); % Take first if multiple matches
                else
                    contrast_p = NaN;
                end
                
                if any(drug_idx)
                    drug_p = anova_table.pValue(drug_idx);
                    drug_p = drug_p(1);
                else
                    drug_p = NaN;
                end
                
                if any(interaction_idx)
                    interaction_p = anova_table.pValue(interaction_idx);
                    interaction_p = interaction_p(1);
                else
                    interaction_p = NaN;
                end
                
                % Post-hoc tests if requested
                if do_posthoc && any([contrast_p, drug_p, interaction_p] < alpha)
                    try
                        anova_results.(field_name).posthoc_contrast = multcompare(rm, 'Contrast');
                        anova_results.(field_name).posthoc_drug = multcompare(rm, 'Drug');
                    catch
                        if display_results
                            fprintf('Post-hoc tests failed for %s\n', field_name);
                        end
                    end
                end
                
                % Add to summary
                summary_data{row_idx, 1} = pop_name;
                summary_data{row_idx, 2} = sizes(iSize);
                summary_data{row_idx, 3} = contrast_p;
                summary_data{row_idx, 4} = drug_p;
                summary_data{row_idx, 5} = interaction_p;
                summary_data{row_idx, 6} = NaN; % No size effect for individual sizes
                summary_data{row_idx, 7} = size(cell_data, 1);
                row_idx = row_idx + 1;
                
                if display_results
                    fprintf('\n--- %s, Size %d ---\n', pop_name, sizes(iSize));
                    fprintf('N cells: %d\n', size(cell_data, 1));
                    fprintf('Contrast effect: p = %.4f %s\n', contrast_p, getSignificanceString(contrast_p, alpha));
                    fprintf('Drug effect: p = %.4f %s\n', drug_p, getSignificanceString(drug_p, alpha));
                    fprintf('Contrast × Drug: p = %.4f %s\n', interaction_p, getSignificanceString(interaction_p, alpha));
                end
                
            catch ME
                if display_results
                    fprintf('ANOVA failed for %s, Size %d: %s\n', pop_name, sizes(iSize), ME.message);
                end
            end
        end
    end
end

% Create summary table
if ~isempty(summary_data)
    stats_table = cell2table(summary_data, ...
        'VariableNames', {'Population', 'Size', 'Contrast_p', 'Drug_p', 'ContrastDrug_p', 'Size_p', 'N_cells'});
    
    if display_results
        fprintf('\n=== SUMMARY TABLE ===\n');
        disp(stats_table);
        
        % Count significant effects
        sig_contrast = sum(stats_table.Contrast_p < alpha);
        sig_drug = sum(stats_table.Drug_p < alpha);
        sig_interaction = sum(stats_table.ContrastDrug_p < alpha);
        
        fprintf('\nSignificant effects (p < %.3f):\n', alpha);
        fprintf('Contrast effects: %d/%d\n', sig_contrast, height(stats_table));
        fprintf('Drug effects: %d/%d\n', sig_drug, height(stats_table));
        fprintf('Contrast × Drug interactions: %d/%d\n', sig_interaction, height(stats_table));
        
        if include_size_effect && nSizes > 1
            sig_size = sum(~isnan(stats_table.Size_p) & stats_table.Size_p < alpha);
            fprintf('Size effects: %d/%d\n', sig_size, sum(~isnan(stats_table.Size_p)));
        end
    end
else
    stats_table = table();
end

% Save results if requested
if ~isempty(save_path)
    save(save_path, 'anova_results', 'stats_table');
    if display_results
        fprintf('\nResults saved to: %s\n', save_path);
    end
end

end

% Helper function to get significance string
function sig_str = getSignificanceString(p_val, alpha)
    if p_val < alpha/1000
        sig_str = '***';
    elseif p_val < alpha/100
        sig_str = '**';
    elseif p_val < alpha
        sig_str = '*';
    else
        sig_str = '';
    end
end