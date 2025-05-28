%% Figure 1
%% Figure 1D - example cell traces and- normalized amplitude plot
load('Figure1D_data.mat')
figure;
subplot(1,3,1); hold on;
plot(timeVector,examplePyr(:,1),'color','k');
plot(timeVector,examplePyr(:,2),'color',[0 0 1]);
fix_axes(gcf,10,'ms','pA'); xlim([-5 20]); ylim([-400 5]); axis square;

subplot(1,3,2); hold on;
plot(timeVector,exampleSST(:,1),'color','k');
plot(timeVector,exampleSST(:,2),'color',[0 0 1]);
fix_axes(gcf,10,'ms','pA'); xlim([-5 20]);ylim([-400 5]); axis square;

subplot(1,3,3); hold on;
plot(1:3,abs([IN_amp_baseline;IN_amp_YM90K;IN_amp_NBQX]./IN_amp_baseline)','color',[1 0.5 0.5]);
plot(1:3,abs([pyr_amp_baseline;pyr_amp_YM90K;pyr_amp_NBQX]./pyr_amp_baseline)','color',[0.6 0.6 0.6]);
fast_errbar_plotting(1:3,abs([IN_amp_baseline;IN_amp_YM90K;IN_amp_NBQX]./IN_amp_baseline),2,'color',[1 0.5 0.5]*0.8);
fast_errbar_plotting(1:3,abs([pyr_amp_baseline;pyr_amp_YM90K;pyr_amp_NBQX]./pyr_amp_baseline),2,'color',[0 0 0]);

xticks(1:3);xticklabels({'Ctrl','+DART','+NBQX'})
ylim([0 1.5]); xlim([0.4 3.5]); fix_axes(gcf,8,'Condition','Norm. Amplitude'); axis square;


%% Figure 1D - ANOVA and pairwise t-tests
load('Figure1_data.mat')
IN_combined=[IN_amp_baseline;IN_amp_YM90K;IN_amp_NBQX]';
pyr_combined=[pyr_amp_baseline;pyr_amp_YM90K;pyr_amp_NBQX]';
%% Mixed ANOVA for drug condition (within) x cell type (between)

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