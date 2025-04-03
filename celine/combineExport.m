% Assuming your four lists are already defined:
% red_avg_stat, red_se_stat, red_avg_loc, red_se_loc
% Each is a vector with 4 elements

% Create a table with 4 rows and 4 columns
combined_table = table(red_avg_stat', red_se_stat', red_avg_loc', red_se_loc');

% Set meaningful variable names for the columns
combined_table.Properties.VariableNames = {'AvgStat', 'SEStat', 'AvgLoc', 'SELoc'};

% If you want, you can also add row names
combined_table.Properties.RowNames = {'12.5%', '25%', '50%', '100%'};

% Save the table as a CSV file
writetable(combined_table, 'PV_data_combined.csv', 'WriteRowNames', true);

disp('4x4 table saved as PV_data_combined.csv');