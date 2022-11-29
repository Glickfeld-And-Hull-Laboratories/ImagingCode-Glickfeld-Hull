function data_struct = TC_stats(trial_averaged_data, base_win, resp_win)        
    %TC mean, std, sem, cell count and more stats for a time X cell response matrix
    
    data_mean = mean(trial_averaged_data,2);
    data_std = std(trial_averaged_data, 0, 2);
    data_count = size(trial_averaged_data,2);
    data_sem = data_std / sqrt(data_count);
    
    data_struct.count_cells = data_count;
    data_struct.mean_TC = data_mean;
    data_struct.std_TC = data_std;
    data_struct.sem_TC = data_sem;
    
    
    %base_win
    base_win_df_by_cell = trial_averaged_data(base_win, :); % Responses in baseline window by cell
    base_win_df_by_cell_mean = mean(base_win_df_by_cell, 1); % average response in window by cell (save this and below)
   
    data_struct.base_win_df_by_cell_mean = base_win_df_by_cell_mean;
    
    %resp_win
    resp_win_df_by_cell = trial_averaged_data(resp_win, :); % Responses in response window by cell
    resp_win_df_by_cell_mean = mean(resp_win_df_by_cell, 1); % average response in window by cell (save this and below)

    data_struct.resp_win_df_by_cell_mean = resp_win_df_by_cell_mean;
    
    % delta f by cell - cell mean of resp - mean of base win
    delta_df = resp_win_df_by_cell_mean - base_win_df_by_cell_mean;

        
    %stats
    delta_df_mean = mean(delta_df); % grand average 
    delta_df_std = std(delta_df);
    delta_df_sem = delta_df_std / sqrt(data_count);
    
    data_struct.base_win = base_win;
    data_struct.resp_win = resp_win;


    data_struct.delta_df_by_cell = delta_df;
    data_struct.mean_delta_df = delta_df_mean;
    data_struct.std_delta_df = delta_df_std;
    data_struct.sem_delta_df = delta_df_sem;      


    %delta f for each adaptor

    %Basline/Responsive windows for each adaptor

    if size(trial_averaged_data,1) == 100

        base_win_each_adapt = zeros(4, length(base_win));
        resp_win_each_adapt = zeros(4, length(resp_win));

        for i = 2:4
            base_win_each_adapt (1,:) = base_win;
            base_win_each_adapt(i,:) = base_win_each_adapt(i-1,:) + 11;

            resp_win_each_adapt (1,:) = resp_win;
            resp_win_each_adapt(i,:) = resp_win_each_adapt(i-1,:) + 11;
        end
    
        data_struct.base_win_each_adapt = base_win_each_adapt;
        data_struct.resp_win_each_adapt = resp_win_each_adapt;


        for i = 1:4
            %base_win
            base_win_each_adapt_df_by_cell = trial_averaged_data(base_win_each_adapt(i,:), :); % Responses in baseline window by cell
            base_win_each_adapt_df_by_cell_mean(i,:) = mean(base_win_each_adapt_df_by_cell, 1); % average response in window by cell (save this and below)
                  
            %resp_win
            resp_win_each_adapt_df_by_cell = trial_averaged_data(resp_win_each_adapt(i,:), :); % Responses in response window by cell
            resp_win_each_adapt_df_by_cell_mean(i,:) = mean(resp_win_each_adapt_df_by_cell, 1); % average response in window by cell (save this and below)
        
            
            % delta f by cell - cell mean of resp - mean of base win
            delta_df_each_adapt_by_cell(i,:) = resp_win_each_adapt_df_by_cell_mean(i,:) - base_win_each_adapt_df_by_cell_mean(i,:);

            %stats
            delta_df_each_adapt_mean(i) = mean(delta_df_each_adapt_by_cell(i,:)); % grand average 
            delta_df_each_adapt_std(i) = std(delta_df_each_adapt_by_cell(i,:));
            delta_df_each_adapt_sem(i) = delta_df_each_adapt_std(i) / sqrt(data_count);         


        end

        Response_A = delta_df_each_adapt_by_cell(1,:);
        Response_B = delta_df_each_adapt_by_cell(2,:);
        
        Aix = Response_B./(Response_A+Response_B);
        Aix_mean = mean(Aix);
        Aix_std = std(Aix);
        Aix_sem = Aix_std / sqrt(length(Aix));

        Aix_alt = Response_B./Response_A;
        Aix_alt_mean = mean(Aix_alt);
        Aix_alt_std = std(Aix_alt);
        Aix_alt_sem = Aix_alt_std / sqrt(length(Aix_alt));


        
        data_struct.base_win_each_adapt_df_by_cell_mean = base_win_each_adapt_df_by_cell_mean;
        data_struct.resp_win_each_adapt_df_by_cell_mean = resp_win_each_adapt_df_by_cell_mean;

        data_struct.delta_df_each_adapt_by_cell = delta_df_each_adapt_by_cell;
        data_struct.mean_delta_df_each_adapt = delta_df_each_adapt_mean;
        data_struct.std_delta_df_each_adapt = delta_df_each_adapt_std;
        data_struct.sem_delta_df_each_adapt = delta_df_each_adapt_sem; 

        data_struct.Aix_by_cell = Aix;
        data_struct.Aix_mean = Aix_mean;
        data_struct.Aix_std = Aix_std;
        data_struct.Aix_sem = Aix_sem;

        data_struct.Aix_alt = Aix_alt;
        data_struct.Aix_alt_mean = Aix_alt_mean;
        data_struct.Aix_alt_std = Aix_alt_std;
        data_struct.Aix_alt_sem = Aix_alt_sem;



    end
    
    
end
