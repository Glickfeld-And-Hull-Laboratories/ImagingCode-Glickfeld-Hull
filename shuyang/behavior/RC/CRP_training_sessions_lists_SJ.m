%CRP_training_days_lists

%% list of variables to initialize

avg_licks_post_cue = [];
avg_licks_pre_cue = [];
avg_licks_post_cue_sem = [];
avg_licks_pre_cue_sem = [];

RT_across_days = [];
RT_across_days_sem = [];
RT_across_days_b2 = [];
RT_across_days_sem_b2 = [];
std_of_RT_across_days = [];
std_of_RT_across_days_b2 = [];

TFT_rates = [];
miss_rates = [];
TFT_rates_b2 = [];
miss_rates_b2 = [];

RT_across_sessions = [];
RT_across_sessions_delay = [];
RT_across_sessions_1000ms_delay = [];
RT_across_sessions_delay_b2 = [];
RT_across_sessions_drift = [];
RT_across_sessions_drift_b2 = [];

days_divider_inx = [];
days_divider_inx_delay = [];
days_divider_inx_1000ms_delay = [];
days_divider_inx_drift = [];
days_divider_inx_delay_b2 = [];

non_consecutive_inx = [];
non_consecutive_inx_delay = [];
non_consecutive_inx_delay_b2 = [];
non_consecutive_inx_1000ms_delay = [];
non_consecutive_inx_drift = [];
non_consecutive_inx_drift_b2 = [];

pre_cue_lick_window_avg = [];
pre_cue_lick_rate_sem = [];
iti_lick_window_avg = [];
iti_lick_rate_sem = [];

%sessions = {'201107_img1073','201108_img1073','201109_img1073','201110_img1073','201111_img1073','201112_img1073'};
%sessions = {'201107_img1074','201108_img1074','201109_img1074','201110_img1074','201111_img1074','201112_img1074'};
sessions = {'201112_img1074-1006','201112_img1074-1033'};

