%% Batch processing script for two-channel imaging data
% This script processes multiple images at different depths and wavelengths
% and compiles all red images and green images into separate PDF pages

% Define the mouse ID and date
mouse = 'i3323';
date = '250319';

% Define the image parameters for each acquisition
% Format: {redFolder, greenFolder, depth, redWavelength, greenWavelength}
acquisitions = {
    {'000', '000', '182, 1x, 4%', '800', '800'},
    {'001', '001', '182, 1x, 18%', '980', '980'},
    {'002', '002', '182, 1x, 46%', '1040', '1040'}
    % Add more acquisitions as needed
};
% Base path for outputs (adjust as needed)
if strcmp(computer, 'GLNXA64')
    out_base = '/home/cc735@dhe.duke.edu/GlickfeldLabShare/All_Staff/home/ACh/Analysis/2p_analysis/SST_YM90K';
else
    out_base = 'Z:\home\ACh\Analysis\2p_analysis';
end

% Create arrays to store all processed images
red_images = cell(length(acquisitions), 1);
green_images = cell(length(acquisitions), 1);
depth_labels = cell(length(acquisitions), 1);
red_wavelength_labels= cell(length(acquisitions), 1);
green_wavelength_labels= cell(length(acquisitions), 1);

% Process all acquisitions and store the results
for i = 1:length(acquisitions)
    params = acquisitions{i};
    redFolder = params{1};
    greenFolder = params{2};
    depth = params{3};
    redWavelength = params{4};
    greenWavelength = params{5};
    
    fprintf('Processing acquisition %d of %d: depth %s\n', i, length(acquisitions), depth);
    
    % Call the processing function
    [red_img, green_img, ~] = processChannelImages(mouse, date, redFolder, greenFolder, depth, redWavelength, greenWavelength);
    
    % Parse depth to remove magnification
    depth_str = params{3};
    depth_val = strtok(depth_str, ',');
    
    % Store the results
    red_images{i} = red_img;
    green_images{i} = green_img;
    depth_labels{i} = depth_val;
    red_wavelength_labels{i}=redWavelength;
    green_wavelength_labels{i}=greenWavelength;
end

% Create a composite PDF with all red images
figure('Position', [100, 100, 1200, 800]);

% Calculate optimal grid dimensions based on number of acquisitions
num_acquisitions = length(acquisitions);
grid_rows = floor(sqrt(num_acquisitions));
grid_cols = ceil(num_acquisitions/grid_rows);

% Create tiled layout with calculated dimensions
t = tiledlayout(grid_rows, grid_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, ['Red channel images - ' mouse ' - ' date], 'FontSize', 14);

for i = 1:length(red_images)
    nexttile;
    imagesc(red_images{i});
    colormap gray;
    title([depth_labels{i} 'm at ' red_wavelength_labels{i} 'nm']);
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    axis tight;
end

% Save the composite red image PDF
out_path = fullfile(out_base, mouse, date, 'batch_output');
if ~exist(out_path, 'dir')
    mkdir(out_path);
end
print(fullfile(out_path, [date '_' mouse '_all_red_images.pdf']), '-dpdf', '-bestfit');

% Create a composite PDF with all green images
figure('Position', [100, 100, 1200, 800]);

% Use the same grid dimensions as calculated for red images
t = tiledlayout(grid_rows, grid_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, ['Green channel images - ' mouse ' - ' date], 'FontSize', 14);

for i = 1:length(green_images)
    nexttile;
    imagesc(green_images{i});
    colormap gray;
    title([depth_labels{i} 'm at ' green_wavelength_labels{i} 'nm']);
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    axis tight;
end

% Save the composite green image PDF
print(fullfile(out_path, [date '_' mouse '_all_green_images.pdf']), '-dpdf', '-bestfit');

% Create a combined document using the export_fig package (if available)
try
    % Check if export_fig is available
    if exist('append_pdfs', 'file')
        red_pdf = fullfile(out_path, [date '_' mouse '_all_red_images.pdf']);
        green_pdf = fullfile(out_path, [date '_' mouse '_all_green_images.pdf']);
        combined_pdf = fullfile(out_path, [date '_' mouse '_all_channels.pdf']);
        append_pdfs(combined_pdf, red_pdf, green_pdf);
        fprintf('Combined PDF created at: %s\n', combined_pdf);
    else
        fprintf('export_fig package not found. Combined PDF not created.\n');
        fprintf('To create a combined PDF, install export_fig from MATLAB File Exchange.\n');
    end
catch
    fprintf('Error creating combined PDF. Individual PDFs are still available.\n');
end

fprintf('Batch processing complete. Results saved to: %s\n', out_path);