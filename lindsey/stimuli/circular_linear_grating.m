function grating_matrix = circular_linear_grating(vertical_size_px, grating_diameter_deg, spatial_freq_cpd, orientation, michelson_contrast, noise_contrast,circle_mask)
    % grating_diameter_deg: diameter of the grating (in degrees of visual angle)
    % spatial_freq_cpd: spatial frequency (cycles per degree)
    % orientation: orientation angle in degrees
    % michelson_contrast: Michelson contrast of the grating (range: 0 to 1)
    % noise_contrast: contrast of the Gaussian white noise (range: 0 to 1)
    % circle_mask: add mask (0 or 1)
    
    % Define constants
    vertical_extent_deg = 80;             % Vertical extent in degrees of visual angle
    aspect_ratio = 16/9;                  % Aspect ratio for standard LCD monitor (16:9)
    horizontal_size_px = round(vertical_size_px * aspect_ratio);  % Adjust horizontal size based on aspect ratio
    mean_luminance = 0.5;                 % Mean luminance for the background (gray)
    
    % Convert degrees to pixels
    deg_to_px = vertical_size_px / vertical_extent_deg;  % Conversion factor: pixels per degree
    grating_diameter_px = round(grating_diameter_deg * deg_to_px);  % Grating diameter in pixels
    
    % Convert orientation to radians
    orientation_rad = deg2rad(orientation);

    % Create a meshgrid that spans the diameter of the grating in pixels
    [X, Y] = meshgrid(linspace(-grating_diameter_px/2, grating_diameter_px/2, grating_diameter_px), ...
                      linspace(-grating_diameter_px/2, grating_diameter_px/2, grating_diameter_px));

    % Rotate the grid based on the orientation
    X_rot = X * cos(orientation_rad) + Y * sin(orientation_rad);
    
    % Define the linear grating using spatial frequency (in cycles per degree)
    % Convert spatial frequency from cycles per degree to cycles per pixel (cpd to cpp)
    spatial_freq_cpp = spatial_freq_cpd / deg_to_px;
    grating = cos(2 * pi * spatial_freq_cpp * X_rot);
    
    % Apply Michelson contrast scaling
    grating = grating * michelson_contrast;

    % Create a circular mask to bound the grating within the defined diameter
    if circle_mask
        r = sqrt(X.^2 + Y.^2);
        mask = r <= (grating_diameter_px / 2);
    
        % Apply the circular mask to the grating
        grating = grating .* mask;
    end

    % Initialize the output matrix to mean luminance (gray)
    grating_matrix = mean_luminance * ones(vertical_size_px, horizontal_size_px);

    % Calculate the center of the matrix to position the grating
    center_x = round(horizontal_size_px / 2);  % Horizontal center based on new matrix size
    center_y = round(vertical_size_px / 2);    % Vertical center

    % Calculate the bounding box for the grating within the output matrix
    half_grating = round(grating_diameter_px / 2);

    % Determine the valid row and column ranges within the grating and matrix
    row_start = max(1, center_y - half_grating + 1);
    row_end = min(vertical_size_px, center_y + half_grating);
    col_start = max(1, center_x - half_grating + 1);
    col_end = min(horizontal_size_px, center_x + half_grating);

    % Adjust grating ranges to match the valid matrix region
    grating_row_start = max(1, half_grating - (center_y - row_start));
    grating_row_end = grating_row_start + (row_end - row_start);
    grating_col_start = max(1, half_grating - (center_x - col_start));
    grating_col_end = grating_col_start + (col_end - col_start);

    % Embed the grating into the central part of the matrix, clipping if necessary
    grating_matrix(row_start:row_end, col_start:col_end) = grating_matrix(row_start:row_end, col_start:col_end) + ...
        grating(grating_row_start:grating_row_end, grating_col_start:grating_col_end) * 0.5;  % Adjust contrast
    
    % Add Gaussian white noise, scaled by noise contrast
    noise_amplitude = michelson_contrast * noise_contrast;
    noise = noise_amplitude * randn(size(grating_matrix));  % Generate Gaussian noise
    grating_matrix = grating_matrix + noise;  % Add noise to the grating

    % Clip values to ensure they are within valid range [0, 1]
    grating_matrix = max(0, min(1, grating_matrix));

    % Display the grating
    figure;
    imshow(grating_matrix, [],'DisplayRange', [0 1], 'InitialMagnification', 'fit');
    colormap(gray);
    axis off;
    title(['Circular Linear Grating: Diameter=' num2str(grating_diameter_deg) '� (degrees), Spatial Frequency=' num2str(spatial_freq_cpd) ' cpd, Orientation=' num2str(orientation) '�, Aspect Ratio=16:9']);
end
