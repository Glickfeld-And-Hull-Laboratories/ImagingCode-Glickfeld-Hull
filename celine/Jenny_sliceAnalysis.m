%load data
load('Z:\home\jen\forCeline\SOM_YM90K_DART_flowIn_amplitudes.mat')
IN_combined=[IN_amp_baseline;IN_amp_YM90K;IN_amp_NBQX]';
pyr_combined=[pyr_amp_baseline;pyr_amp_YM90K;pyr_amp_NBQX]';
%%
% Mixed ANOVA for drug condition (within) x cell type (between)
% Assuming IN_combined and pyr_combined are matrices with dimensions:
% (number of neurons) x (3 drug conditions: ACSF, DART, NBQX)

% First, let's prepare the data for the ANOVA
% We need to create a table with all observations

% Get sample sizes
n_in = size(IN_combined, 1);  % Number of interneurons
n_pyr = size(pyr_combined, 1); % Number of pyramidal cells
n_total = n_in + n_pyr;

% Combine data
all_data = [IN_combined; pyr_combined];

% Create factors
drug_labels = {'ACSF', 'DART', 'NBQX'};
cell_type = [ones(n_in, 1); 2*ones(n_pyr, 1)]; % 1 for interneurons, 2 for pyramidal
subject_id = (1:n_total)';

% Create a table for the mixed ANOVA
anova_table = table(subject_id, cell_type, 'VariableNames', {'Subject', 'CellType'});
anova_table.CellType = categorical(anova_table.CellType, [1 2], {'Interneuron', 'Pyramidal'});

% Add the drug conditions as additional variables
for i = 1:3
    anova_table.(drug_labels{i}) = all_data(:, i);
end

% Create a within-subjects design table
within_design = table(categorical(drug_labels', 'Ordinal', false), 'VariableNames', {'Drug'});

% Setup the repeated measures model
rm = fitrm(anova_table, 'ACSF,DART,NBQX ~ CellType', 'WithinDesign', within_design);

% Run the ANOVA for within-subjects effects (drug and drug*celltype)
ranova_results = ranova(rm);
disp('Within-Subjects ANOVA Results:');
disp(ranova_results);

% Get between-subjects effects (cell type)
between_results = anova(rm);
disp('Between-Subjects ANOVA Results (Cell Type):');
disp(between_results);


% Manually perform pairwise t-tests within each cell type with Bonferroni correction
% Number of comparisons per cell type
num_comparisons = 3;

disp('Pairwise comparisons for drug conditions within each cell type (with Bonferroni correction):');
disp('Interneurons:');

% Interneuron comparisons
fprintf('ACSF vs DART: ');
[~, p_in_acsf_dart] = ttest(IN_combined(:,1), IN_combined(:,2));
p_in_acsf_dart_corr = min(p_in_acsf_dart * num_comparisons, 1); % Apply Bonferroni correction, cap at 1
fprintf('p = %.4f (corrected p = %.4f)\n', p_in_acsf_dart, p_in_acsf_dart_corr);

fprintf('ACSF vs NBQX: ');
[~, p_in_acsf_nbqx] = ttest(IN_combined(:,1), IN_combined(:,3));
p_in_acsf_nbqx_corr = min(p_in_acsf_nbqx * num_comparisons, 1);
fprintf('p = %.4f (corrected p = %.4f)\n', p_in_acsf_nbqx, p_in_acsf_nbqx_corr);

fprintf('DART vs NBQX: ');
[~, p_in_dart_nbqx] = ttest(IN_combined(:,2), IN_combined(:,3));
p_in_dart_nbqx_corr = min(p_in_dart_nbqx * num_comparisons, 1);
fprintf('p = %.4f (corrected p = %.4f)\n', p_in_dart_nbqx, p_in_dart_nbqx_corr);

% Pyramidal comparisons
disp('Pyramidal cells:');
fprintf('ACSF vs DART: ');
[~, p_pyr_acsf_dart] = ttest(pyr_combined(:,1), pyr_combined(:,2));
p_pyr_acsf_dart_corr = min(p_pyr_acsf_dart * num_comparisons, 1);
fprintf('p = %.4f (corrected p = %.4f)\n', p_pyr_acsf_dart, p_pyr_acsf_dart_corr);

fprintf('ACSF vs NBQX: ');
[~, p_pyr_acsf_nbqx] = ttest(pyr_combined(:,1), pyr_combined(:,3));
p_pyr_acsf_nbqx_corr = min(p_pyr_acsf_nbqx * num_comparisons, 1);
fprintf('p = %.4f (corrected p = %.4f)\n', p_pyr_acsf_nbqx, p_pyr_acsf_nbqx_corr);

fprintf('DART vs NBQX: ');
[~, p_pyr_dart_nbqx] = ttest(pyr_combined(:,2), pyr_combined(:,3));
p_pyr_dart_nbqx_corr = min(p_pyr_dart_nbqx * num_comparisons, 1);
fprintf('p = %.4f (corrected p = %.4f)\n', p_pyr_dart_nbqx, p_pyr_dart_nbqx_corr);