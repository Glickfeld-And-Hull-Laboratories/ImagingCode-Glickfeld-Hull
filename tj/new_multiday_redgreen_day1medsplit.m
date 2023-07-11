day1_green_grand_med = median(green_max_all(1,:));

day1_green_bot_50_d_score_prefori_d1_d2 = green_d_score_prefori_d1_d2(green_max_all(1,:) < day1_green_grand_med);
day1_green_top_50_d_score_prefori_d1_d2 = green_d_score_prefori_d1_d2(green_max_all(1,:) > day1_green_grand_med);

day1_green_bot_50_d_score_prefori_d1_d3 = green_d_score_prefori_d1_d3(green_max_all(1,:) < day1_green_grand_med);
day1_green_top_50_d_score_prefori_d1_d3 = green_d_score_prefori_d1_d3(green_max_all(1,:) > day1_green_grand_med);

day1_red_grand_med = median(day1_red_grand_max);

day1_red_bot_50_d_score_prefori_d1_d2 = red_d_score_prefori_d1_d2(day1_red_grand_max < day1_red_grand_med);
day1_red_top_50_d_score_prefori_d1_d2 = red_d_score_prefori_d1_d2(day1_red_grand_max > day1_red_grand_med);

day1_red_bot_50_d_score_prefori_d1_d3 = red_d_score_prefori_d1_d3(day1_red_grand_max < day1_red_grand_med);
day1_red_top_50_d_score_prefori_d1_d3 = red_d_score_prefori_d1_d3(day1_red_grand_max > day1_red_grand_med);

if str2double(expt(d1).img_day) == 1

    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_day1_median_sep_d_scores.mat']), 'day1_green_bot_50_d_score_prefori_d1_d2', 'day1_green_top_50_d_score_prefori_d1_d2', ...
        'day1_green_bot_50_d_score_prefori_d1_d3', 'day1_green_top_50_d_score_prefori_d1_d3', 'day1_red_bot_50_d_score_prefori_d1_d2', 'day1_red_top_50_d_score_prefori_d1_d2', ...
        'day1_red_bot_50_d_score_prefori_d1_d3', 'day1_red_top_50_d_score_prefori_d1_d3')
else
    save(fullfile(newfnout, [mouse, '_', img_area, '_', img_layer, '_day1_median_sep_d_scores_2.mat']), 'day1_green_bot_50_d_score_prefori_d1_d2', 'day1_green_top_50_d_score_prefori_d1_d2', ...
        'day1_green_bot_50_d_score_prefori_d1_d3', 'day1_green_top_50_d_score_prefori_d1_d3', 'day1_red_bot_50_d_score_prefori_d1_d2', 'day1_red_top_50_d_score_prefori_d1_d2', ...
        'day1_red_bot_50_d_score_prefori_d1_d3', 'day1_red_top_50_d_score_prefori_d1_d3')
end

%%

