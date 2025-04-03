function exportLMEResultsTable(lme_sst, lme_pyr, custom_filename)
    % Get the model names from the input arguments
    sst_model_name = inputname(1);
    pyr_model_name = inputname(2);
    
    % If no names are available (e.g., if models were passed as outputs from other functions),
    % use generic model names
    if isempty(sst_model_name)
        sst_model_name = 'lme_sst';
    end
    if isempty(pyr_model_name)
        pyr_model_name = 'lme_pyr';
    end
    
    % Generate filename based on model names
    if nargin < 3 || isempty(custom_filename)
        % Extract meaningful parts from model names, e.g., "lme_sst_stat" -> "stat"
        sst_suffix = regexp(sst_model_name, '_([^_]+)$', 'tokens');
        if ~isempty(sst_suffix)
            model_type = sst_suffix{1}{1};
        else
            model_type = 'model';
        end
        
        filename = [model_type '_results.csv'];
    else
        filename = custom_filename;
    end
    
    % Extract fixed effects from both models
    fe_sst = lme_sst.Coefficients;
    fe_pyr = lme_pyr.Coefficients;
    
    % Create tables for each model
    sst_table = table(fe_sst.Name, fe_sst.Estimate, fe_sst.SE, fe_sst.pValue, ...
                     'VariableNames', {'Parameter', 'Estimate_SST', 'SE_SST', 'pValue_SST'});
                 
    pyr_table = table(fe_pyr.Name, fe_pyr.Estimate, fe_pyr.SE, fe_pyr.pValue, ...
                     'VariableNames', {'Parameter', 'Estimate_PYR', 'SE_PYR', 'pValue_PYR'});
    
    % Join the tables
    params = unique([fe_sst.Name; fe_pyr.Name]);
    combined_table = table(params, 'VariableNames', {'Parameter'});
    
    % Fill in SST values
    combined_table.Estimate_SST = NaN(height(combined_table), 1);
    combined_table.SE_SST = NaN(height(combined_table), 1);
    combined_table.pValue_SST = NaN(height(combined_table), 1);
    
    % Fill in PYR values
    combined_table.Estimate_PYR = NaN(height(combined_table), 1);
    combined_table.SE_PYR = NaN(height(combined_table), 1);
    combined_table.pValue_PYR = NaN(height(combined_table), 1);
    
    % Match parameters and fill in values
    for i = 1:height(combined_table)
        param = combined_table.Parameter{i};
        
        % SST model
        idx_sst = find(strcmp(fe_sst.Name, param));
        if ~isempty(idx_sst)
            combined_table.Estimate_SST(i) = fe_sst.Estimate(idx_sst);
            combined_table.SE_SST(i) = fe_sst.SE(idx_sst);
            combined_table.pValue_SST(i) = fe_sst.pValue(idx_sst);
        end
        
        % PYR model
        idx_pyr = find(strcmp(fe_pyr.Name, param));
        if ~isempty(idx_pyr)
            combined_table.Estimate_PYR(i) = fe_pyr.Estimate(idx_pyr);
            combined_table.SE_PYR(i) = fe_pyr.SE(idx_pyr);
            combined_table.pValue_PYR(i) = fe_pyr.pValue(idx_pyr);
        end
    end
    
    % Add formatted p-values with significance indicators
    combined_table.pValue_SST_Fmt = cell(height(combined_table), 1);
    combined_table.pValue_PYR_Fmt = cell(height(combined_table), 1);
    
    for j = 1:height(combined_table)
        % SST p-values
        if ~isnan(combined_table.pValue_SST(j))
            p = combined_table.pValue_SST(j);
            if p < 0.001
                combined_table.pValue_SST_Fmt{j} = sprintf('%.4f***', p);
            elseif p < 0.01
                combined_table.pValue_SST_Fmt{j} = sprintf('%.4f**', p);
            elseif p < 0.05
                combined_table.pValue_SST_Fmt{j} = sprintf('%.4f*', p);
            else
                combined_table.pValue_SST_Fmt{j} = sprintf('%.4f', p);
            end
        else
            combined_table.pValue_SST_Fmt{j} = '';
        end
        
        % PYR p-values
        if ~isnan(combined_table.pValue_PYR(j))
            p = combined_table.pValue_PYR(j);
            if p < 0.001
                combined_table.pValue_PYR_Fmt{j} = sprintf('%.4f***', p);
            elseif p < 0.01
                combined_table.pValue_PYR_Fmt{j} = sprintf('%.4f**', p);
            elseif p < 0.05
                combined_table.pValue_PYR_Fmt{j} = sprintf('%.4f*', p);
            else
                combined_table.pValue_PYR_Fmt{j} = sprintf('%.4f', p);
            end
        else
            combined_table.pValue_PYR_Fmt{j} = '';
        end
    end
    
    % Create estimate columns with standard errors
    combined_table.Estimate_SE_SST = cell(height(combined_table), 1);
    combined_table.Estimate_SE_PYR = cell(height(combined_table), 1);
    
    for k = 1:height(combined_table)
        if ~isnan(combined_table.Estimate_SST(k))
            combined_table.Estimate_SE_SST{k} = sprintf('%.4f ± %.4f', combined_table.Estimate_SST(k), combined_table.SE_SST(k));
        else
            combined_table.Estimate_SE_SST{k} = '';
        end
        
        if ~isnan(combined_table.Estimate_PYR(k))
            combined_table.Estimate_SE_PYR{k} = sprintf('%.4f ± %.4f', combined_table.Estimate_PYR(k), combined_table.SE_PYR(k));
        else
            combined_table.Estimate_SE_PYR{k} = '';
        end
    end
    
    % Create export table with formatted columns
    export_table = table(combined_table.Parameter, ...
                         combined_table.Estimate_SE_SST, ...
                         combined_table.pValue_SST_Fmt, ...
                         combined_table.Estimate_SE_PYR, ...
                         combined_table.pValue_PYR_Fmt, ...
                         'VariableNames', {'Parameter', 'Estimate_SST', 'pValue_SST', 'Estimate_PYR', 'pValue_PYR'});
    
    % Add model information
    model_info = sprintf('%s vs %s | SST: AIC=%.2f, BIC=%.2f | PYR: AIC=%.2f, BIC=%.2f', ...
                         sst_model_name, pyr_model_name, ...
                         lme_sst.ModelCriterion.AIC, lme_sst.ModelCriterion.BIC, ...
                         lme_pyr.ModelCriterion.AIC, lme_pyr.ModelCriterion.BIC);
    
    % Display the table and model info
    disp('============= Combined Mixed Effects Model Results =============');
    disp(export_table);
    disp('Significance codes: *** p<0.001, ** p<0.01, * p<0.05');
    disp(model_info);
    
    % Export to CSV
    writetable(export_table, filename);
    
    % Add model info to the top of the CSV
    % Read the CSV content
    fid = fopen(filename, 'r');
    content = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    content = content{1};
    
    % Write the content back with model info at the top
    fid = fopen(filename, 'w');
    fprintf(fid, '%s\n', model_info);
    fprintf(fid, 'Significance codes: *** p<0.001, ** p<0.01, * p<0.05\n');
    fprintf(fid, '\n');
    for i = 1:length(content)
        fprintf(fid, '%s\n', content{i});
    end
    fclose(fid);
    
    % Also create a second CSV with the raw values for further analysis
    raw_export_table = table(combined_table.Parameter, ...
                           combined_table.Estimate_SST, combined_table.SE_SST, combined_table.pValue_SST, ...
                           combined_table.Estimate_PYR, combined_table.SE_PYR, combined_table.pValue_PYR, ...
                           'VariableNames', {'Parameter', 'Estimate_SST', 'SE_SST', 'pValue_SST', ...
                                           'Estimate_PYR', 'SE_PYR', 'pValue_PYR'});
    
    [filepath, name, ext] = fileparts(filename);
    raw_filename = fullfile(filepath, [name '_raw' ext]);
    writetable(raw_export_table, raw_filename);
    
    fprintf('Results exported to:\n  %s (formatted)\n  %s (raw values)\n', filename, raw_filename);
end