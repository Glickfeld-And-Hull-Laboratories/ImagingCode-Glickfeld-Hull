function analyzeResponsiveCells(respToLarge, respToSmall, respToLarge_red, respToLarge_green, ...
                                red_concat, green_concat, red_ind_concat, green_ind_concat, ...
                                mouseInds, mouseNames)
% ANALYZERESPONSIVECELLS Analyze and visualize distribution of responsive cells by mouse
%
% Usage:
%   analyzeResponsiveCells(respToLarge, respToSmall, respToLarge_red, respToLarge_green, ...
%                          red_concat, green_concat, red_ind_concat, green_ind_concat, ...
%                          mouseInds, mouseNames)
%
% Inputs:
%   respToLarge - Logical array indicating cells responsive to large stimuli
%   respToSmall - Logical array indicating cells responsive to small stimuli  
%   respToLarge_red - Indices of red cells responsive to large stimuli
%   respToLarge_green - Indices of green cells responsive to large stimuli
%   red_concat - Logical array indicating red cells
%   green_concat - Logical array indicating green cells
%   red_ind_concat - Indices of red cells
%   green_ind_concat - Indices of green cells
%   mouseInds - Cell array with neuron indices for each mouse
%   mouseNames - Cell array with mouse ID strings

% Find indices for each response category
respToSmall_red = find(respToSmall.*red_concat');
respToSmall_green = find(respToSmall.*green_concat');

both_indices = find(respToLarge & respToSmall);
respToBoth_red = intersect(both_indices, red_ind_concat);
respToBoth_green = intersect(both_indices, green_ind_concat);

% Calculate counts for each response category
largeOnly = sum(respToLarge & ~respToSmall);
smallOnly = sum(~respToLarge & respToSmall);
both = sum(respToLarge & respToSmall);
neither = sum(~respToLarge & ~respToSmall);

nMice = length(mouseInds);
red_counts_large = zeros(nMice, 1);
green_counts_large = zeros(nMice, 1);
red_counts_small = zeros(nMice, 1);
green_counts_small = zeros(nMice, 1);
red_counts_both = zeros(nMice, 1);
green_counts_both = zeros(nMice, 1);

for i = 1:nMice
    red_counts_large(i) = length(intersect(mouseInds{i}, respToLarge_red));
    green_counts_large(i) = length(intersect(mouseInds{i}, respToLarge_green));
    red_counts_small(i) = length(intersect(mouseInds{i}, respToSmall_red));
    green_counts_small(i) = length(intersect(mouseInds{i}, respToSmall_green));
    red_counts_both(i) = length(intersect(mouseInds{i}, respToBoth_red));
    green_counts_both(i) = length(intersect(mouseInds{i}, respToBoth_green));
end

mouse_labels = mouseNames(:);

% Create table for large-responsive cells
total_counts_large = red_counts_large + green_counts_large;
results_table_large = table(mouse_labels, red_counts_large, green_counts_large, total_counts_large, ...
    'VariableNames', {'Mouse', 'Red_Large_Responsive', 'Green_Large_Responsive', 'Total'});
disp('Distribution of Large-Responsive Cells by Mouse:');
disp(results_table_large);
fprintf('Total red large-responsive: %d\n', sum(red_counts_large));
fprintf('Total green large-responsive: %d\n', sum(green_counts_large));

% Create table for small-responsive cells
total_counts_small = red_counts_small + green_counts_small;
results_table_small = table(mouse_labels, red_counts_small, green_counts_small, total_counts_small, ...
    'VariableNames', {'Mouse', 'Red_Small_Responsive', 'Green_Small_Responsive', 'Total'});
disp('Distribution of Small-Responsive Cells by Mouse:');
disp(results_table_small);
fprintf('Total red small-responsive: %d\n', sum(red_counts_small));
fprintf('Total green small-responsive: %d\n', sum(green_counts_small));

% Create table for "both" cells
total_counts_both = red_counts_both + green_counts_both;
results_table_both = table(mouse_labels, red_counts_both, green_counts_both, total_counts_both, ...
    'VariableNames', {'Mouse', 'Red_Both_Responsive', 'Green_Both_Responsive', 'Total'});
disp('Distribution of Both-Responsive Cells by Mouse:');
disp(results_table_both);
fprintf('Total red both-responsive: %d\n', sum(red_counts_both));
fprintf('Total green both-responsive: %d\n', sum(green_counts_both));

% Create pie charts in 3x2 layout
figure;

subplot(3,2,1)
pie(red_counts_large);
title({'Red Large-Responsive Cells', sprintf('(n=%d)', sum(red_counts_large))});
set(gca, 'TickDir', 'out');
grid off;
box off;

subplot(3,2,2)
pie(green_counts_large);
title({'Green Large-Responsive Cells', sprintf('(n=%d)', sum(green_counts_large))});
set(gca, 'TickDir', 'out');
grid off;
box off;

subplot(3,2,3)
pie(red_counts_small);
title({'Red Small-Responsive Cells', sprintf('(n=%d)', sum(red_counts_small))});
set(gca, 'TickDir', 'out');
grid off;
box off;

subplot(3,2,4)
pie(green_counts_small);
title({'Green Small-Responsive Cells', sprintf('(n=%d)', sum(green_counts_small))});
set(gca, 'TickDir', 'out');
grid off;
box off;

subplot(3,2,5)
pie(red_counts_both);
title({'Red Both-Responsive Cells', sprintf('(n=%d)', sum(red_counts_both))});
set(gca, 'TickDir', 'out');
grid off;
box off;

subplot(3,2,6)
pie(green_counts_both);
title({'Green Both-Responsive Cells', sprintf('(n=%d)', sum(green_counts_both))});
set(gca, 'TickDir', 'out');
grid off;
box off;

x0=5;
y0=5;
width=9;
height=6;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')
end