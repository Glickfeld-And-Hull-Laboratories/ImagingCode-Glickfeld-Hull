function [gOSI_output] = gOSI(input_df_matrix,oris)
%gOSI calculate global osi of cells
%   Based on code from Lan Lou, 2023, which itself is based on https://openreview.net/pdf?id=Bkx_Dj09tQ
%   Input should be a matrix of [df/f] responses in the shape nCells x
%   nDirections/orientations
gOSI_output=nan(size(1,length(input_df_matrix)));
    for iCell = 1:length(input_df_matrix)
        
        df_cell = input_df_matrix(iCell, :);
        
        tuning = df_cell - min(df_cell); % ensure all values are non-negative
        
        theta_arr = oris; % according to formula: unit is deg, not rad
        sin_arr = sin(2 * deg2rad(theta_arr));
        cos_arr = cos(2 * deg2rad(theta_arr));
        
        gOSI= sqrt((sum(tuning .* sin_arr))^2 + (sum(tuning .* cos_arr))^2) / sum(tuning);
        
        gOSI_output(iCell)=gOSI;
    end

end