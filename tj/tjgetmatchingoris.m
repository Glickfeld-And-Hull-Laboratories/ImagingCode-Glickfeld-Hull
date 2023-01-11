% function [prefori_d1_d2_match_tune, prefori_d2_match_tune, d_scores_all, k_d1_d2_match, k_d2_match, dscore_k_d1_d2_match, max_d1_d2_match, max_d2_match, dscore_max_d1_d2_match] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max)
function [a, b, c, d, e, f, g, h, i] = tjgetmatchingoris(d1_ori, d2_ori, d2_matches, d1_k_max, d2_k_max)

%find tuned cells for each day and cells that are matched across days adn
%tuned on at least one
tuned_d1 = d1_ori.ind_theta90;
tuned_d2 = d2_ori.ind_theta90;

match_d2 = find([d2_matches.cellImageAlign.pass]); 

matched_d1_tuned_d2 = intersect(match_d2, tuned_d1);
matched_d1_d2_tuned = intersect(match_d2, tuned_d2);
matched_d1_d2_all_tuned = unique([matched_d1_tuned_d2 matched_d1_d2_tuned]);


prefori_d1_d2_match = d1_ori.prefOri(1,match_d2);
prefori_d2_match = d2_ori.prefOri(1,match_d2);


%find pref oris of green matched cells for d1 and d2
prefori_d1_d2_match_tune = d1_ori.prefOri(1,matched_d1_d2_all_tuned);
prefori_d2_match_tune = d2_ori.prefOri(1,matched_d1_d2_all_tuned);


%%
%diff scores for pref ori

%d1_2
dscore_prefori_d1_d2_match = double(prefori_d1_d2_match_tune>90);
dscore_prefori_d1_d2_match(dscore_prefori_d1_d2_match>0) = 180;
dscore_prefori_d1_d2_match = abs(dscore_prefori_d1_d2_match-prefori_d1_d2_match_tune);
%d2
dscore_prefori_d2_match = double(prefori_d2_match_tune>90);
dscore_prefori_d2_match(dscore_prefori_d2_match>0) = 180;
dscore_prefori_d2_match = abs(dscore_prefori_d2_match-prefori_d2_match_tune);

d_score_prefori_d1_d2 = abs(dscore_prefori_d1_d2_match-dscore_prefori_d2_match);

d_scores_all = [d_score_prefori_d1_d2];

%% K VALUES

%find ks of matched cells for d1 and d2
k_d1_d2_match = d1_k_max.k1(1,matched_d1_d2_all_tuned);
k_d2_match = d2_k_max.k1(1,matched_d1_d2_all_tuned);

%d1_2
dscore_k_d1_d2_match = abs(k_d1_d2_match-k_d2_match);


%% MAX VALUES


%find max of matched cells for d1 and d2
max_d1_d2_match = d1_k_max.max_dfof(1,matched_d1_d2_all_tuned);
max_d2_match = d2_k_max.max_dfof(1,matched_d1_d2_all_tuned);



%d1_2
dscore_max_d1_d2_match = abs(max_d1_d2_match-max_d2_match);


%%
a = prefori_d1_d2_match_tune;
b = prefori_d2_match_tune;
c = d_scores_all;
d = k_d1_d2_match;
e = k_d2_match;
f = dscore_k_d1_d2_match;
g = max_d1_d2_match;
h = max_d2_match;
i = dscore_max_d1_d2_match;
end