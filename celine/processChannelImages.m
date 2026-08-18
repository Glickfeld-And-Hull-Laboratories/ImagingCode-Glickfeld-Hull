function [redImage, greenImage, compositeImage] = processChannelImages(mouse, date, redFolder, greenFolder, depth, redWavelength, greenWavelength)
% PROCESSCHANNELIMAGES Processes two-channel microscopy images and returns processed images
%   [redImage, greenImage, compositeImage] = processChannelImages(mouse, date, redFolder, 
%   greenFolder, depth, redWavelength, greenWavelength)
%
%   Inputs:
%   mouse - Mouse ID (e.g., 'i2184')
%   date - Date of recording (e.g., '241217')
%   redFolder - First three digits of red channel folder (e.g., '001')
%   greenFolder - First three digits of green channel folder (e.g., '000')
%   depth - Imaging depth information (e.g., '133, 2x')
%   redWavelength - (Optional) Wavelength used for red channel (default: '1040')
%   greenWavelength - (Optional) Wavelength used for green channel (default: '920')
%
%   Outputs:
%   redImage - Processed red channel image
%   greenImage - Processed green channel image
%   compositeImage - Composite of red and green channels

    % Set default values for optional parameters
    if nargin < 7 || isempty(greenWavelength)
        greenWavelength = '920';
    end
    
    if nargin < 6 || isempty(redWavelength)
        redWavelength = '1040';
    end

    % Determine base paths based on computer type
    if strcmp(computer, 'GLNXA64')
        base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\Bonnie\Data';
        out_base = '/home/cc735@dhe.duke.edu/GlickfeldLabShare/All_Staff/home/ACh/Analysis/2p_analysis/SST_YM90K';
    % else if strcmp(computer, 'PCWIN64')
    %     base = 'G:\home\ACh\Data\2p_data';
    %     out_base = 'G:\home\ACh\Analysis\2p_analysis';
    else
        base = '\home\ACh\Data\2p_data';
        out_base = '\home\ACh\Analysis\2p_analysis';
    end

    % Create paths
    data_path_red = fullfile(base, mouse, date, redFolder);
    out_path = fullfile(out_base, mouse, date, redFolder);
    
    % Create output directory if it doesn't exist
    if ~exist(out_path, 'dir')
        mkdir(out_path);
    end
    
    % Process red channel
    cd(data_path_red);
    nframes_red = 1000;
    % Hardcoded redrun as '000'
    load([redFolder '_000_000.mat']);
    data = sbxread([redFolder '_000_000'], 0, nframes_red);
    
    data_g = squeeze(data(1,:,:,:));
    data_r = squeeze(data(2,:,:,:));
    
    data_g_avg = mean(data_g(:,:,:), 3);
    data_r_avg = mean(data_r(:,:,:), 3);
    
    % Register red channel
    [out, data_r_reg] = stackRegister(data_r, data_r_avg);
    [outs, data_g_reg] = stackRegister_MA(data_g, [], [], out);
    
    data_g_reg_avg = mean(data_g_reg, 3);
    data_r_reg_avg = mean(data_r_reg, 3);
    
    % Save red channel image
    cd(out_path);
    fig1 = figure; 
    imagesc(data_r_reg_avg);
    title([' ' depth ' red at ' redWavelength]);
    colormap gray;
    %caxis([200 1000])
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    print(fullfile(out_path, [date '_' mouse '_red_FOV.pdf']), '-dpdf', '-bestfit');
    saveas(fig1, 'redFOV.png');
    
    % Process green channel
    nframes_green = 200;
    data_path_green = fullfile(base, mouse, date, greenFolder);
    cd(data_path_green);
    
    % Hardcoded greenrun as '000'
    load([greenFolder '_000_000.mat']);
    data_920 = sbxread([greenFolder '_000_000'], 0, nframes_green);
    
    data_g_920 = squeeze(data_920(1,:,:,:));
    data_g_920_avg = mean(data_g_920(:,:,:), 3);
    
    [out, data_g_reg_920] = stackRegister(data_g_920, data_g_reg_avg);
    regImg = mean(data_g_reg_920, 3);
    
    % Save green channel image
    fig2 = figure;
    imagesc(regImg);
    colormap gray;
    caxis([200 4000]);
    cd(out_path);
    title([' ' depth ' green at ' greenWavelength]);
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    print(fullfile(out_path, [date '_' mouse '_FOV.pdf']), '-dpdf', '-bestfit');
    saveas(fig2, 'greenFOV.png');
    
    % Create and save composite image
    sz = size(regImg);
    rgb = zeros(sz(1), sz(2), 3);
    rgb(:,:,1) = data_r_reg_avg ./ max(data_r_reg_avg(:));
    rgb(:,:,2) = regImg ./ max(regImg(:));
    
    fig3 = figure; 
    image(rgb);
    title([' ' depth]);
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    print(fullfile(out_path, [date '_' mouse 'red_and_green_FOV.pdf']), '-dpdf', '-bestfit');
    saveas(fig3, 'redAndgreenFOV.png');
    
    % Return processed images
    redImage = data_r_reg_avg;
    greenImage = regImg;
    compositeImage = rgb;
    
    % Close figures if no output arguments requested
    if nargout == 0
        close(fig1);
        close(fig2);
        close(fig3);
    end
end