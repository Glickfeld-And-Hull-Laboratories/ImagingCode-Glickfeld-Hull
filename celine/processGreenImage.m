function greenImage = processGreenImage(mouse, date, greenFolder, depth, greenWavelength, varargin)
% PROCESSGREENIMAGE Processes green channel microscopy images and returns processed image
%   greenImage = processGreenImage(mouse, date, greenFolder, depth, greenWavelength, nframes)
%
%   Inputs:
%   mouse - Mouse ID (e.g., 'i2184')
%   date - Date of recording (e.g., '241217')
%   greenFolder - First three digits of green channel folder (e.g., '000')
%   depth - Imaging depth information (e.g., '133, 2x')
%   greenWavelength - Wavelength used for green channel (e.g., '920')
%   nframes - (Optional) Number of frames to process (default: 200)
%
%   Outputs:
%   greenImage - Processed green channel image

    % Determine base paths based on computer type
    if strcmp(computer, 'GLNXA64')
        base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\Bonnie\Data';
        out_base = '/home/cc735@dhe.duke.edu/GlickfeldLabShare/All_Staff/home/ACh/Analysis/2p_analysis/SST_YM90K';
    else
        base = 'Z:\home\ACh\Data\2p_data';
        out_base = 'Z:\home\ACh\Analysis\2p_analysis';
    end

    % Create paths
    data_path_green = fullfile(base, mouse, date, greenFolder);
    out_path = fullfile(out_base, mouse, date, greenFolder);
    
    % Create output directory if it doesn't exist
    if ~exist(out_path, 'dir')
        mkdir(out_path);
    end
    
    % Process green channel
    % Parse input arguments for number of frames
    p = inputParser;
    addOptional(p, 'nframes', 200, @isnumeric);
    parse(p, varargin{:});
    nframes_green = p.Results.nframes;
    
    cd(data_path_green);
    
    % Load green channel data
    load([greenFolder '_000_000.mat']);
    data_green = sbxread([greenFolder '_000_000'], 0, nframes_green);
    
    data_g = squeeze(data_green(1,:,:,:));
    data_g_avg = mean(data_g(:,:,:), 3);
    
    % Register green channel
    [out, data_g_reg] = stackRegister(data_g, data_g_avg);
    greenImage = mean(data_g_reg, 3);
    
    % Save green channel image
    cd(out_path);
    fig = figure;
    imagesc(greenImage);
    colormap gray;
    %caxis([200 4000]);
    title([' ' depth ' green at ' greenWavelength]);
    set(gca, 'TickDir', 'out');
    grid off;
    box off;
    print(fullfile(out_path, [date '_' mouse '_FOV.pdf']), '-dpdf', '-bestfit');
    saveas(fig, 'greenFOV.png');
    
    % Close figure if no output arguments requested
    if nargout == 0
        close(fig);
    end
end