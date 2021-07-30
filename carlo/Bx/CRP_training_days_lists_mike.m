%CRP_training_days_lists

%% list of variables to initialize
initVar = input('Load bx variable blanks? If so, type 1, else type 0: ');
if initVar == 1
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
else
end

%days = {'191023_img1505','191024_img1505','191025_img1505','191026_img1505','191027_img1505','191028_img1505', '191029_img1505','191030_img1505'}; 
%days = {'191023_img1506','191024_img1506','191025_img1506','191026_img1506','191027_img1506','191028_img1506', '191029_img1506','191030_img1506'};
%days = {'191103_img1507','191104_img1507','191105_img1507','191106_img1507','191107_img1507','191108_img1507','191109_img1507','191110_img1507','191111_img1507','191112_img1507','191113_img1507','191114_img1507','191115_img1507','191118_img1507','191119_img1507'};  
%days = {'191102_img1508','191103_img1508','191104_img1508','191105_img1508','191106_img1508','191107_img1508','191108_img1508','191109_img1508','191110_img1508','191111_img1508','191112_img1508','191113_img1508','191114_img1508','191115_img1508','191118_img1508','191119_img1508','191121_img1508'};
%days = {'200106_img1510','200107_img1510','200108_img1510','200109_img1510','200110_img1510','200111_img1510','200112_img1510','200113_img1510','200115_img1510','200116_img1510','200117_img1510','200119_img1510','200120_img1510','200121_img1510','200122_img1510','200123_img1510','200124_img1510','200125_img1510','200126_img1510','200129_img1510','200130_img1510','200131_img1510','200203_img1510'};
%days = {'200103_img1511','200104_img1511','200105_img1511','200106_img1511','200107_img1511','200108_img1511','200109_img1511','200110_img1511','200111_img1511','200112_img1511','200113_img1511','200115_img1511','200116_img1511','200117_img1511','200119_img1511','200120_img1511','200121_img1511','200122_img1511','200123_img1511','200124_img1511','200125_img1511','200126_img1511'};
%days = {'200310_img1512', '200311_img1512', '200312_img1512', '200313_img1512', '200314_img1512', '200316_img1512', '200317_img1512', '200318_img1512', '200319_img1512', '200320_img1512', '200321_img1512', '200323_img1512', '200324_img1512', '200325_img1512', '200326_img1512', '200327_img1512', '200328_img1512', '200329_img1512', '200330_img1512', '200331_img1512'};
%days = {'200312_img1513', '200313_img1513', '200314_img1513', '200316_img1513', '200317_img1513', '200318_img1513', '200319_img1513', '200320_img1513', '200321_img1513', '200323_img1513', '200324_img1513', '200325_img1513', '200326_img1513', '200327_img1513', '200328_img1513', '200329_img1513', '200330_img1513', '200331_img1513', '200401_img1513', '200402_img1513', '200407_img1513', '200408_img1513', '200409_img1513'};
%days = {'200719_img1516', '200720_img1516', '200721_img1516', '200722_img1516', '200723_img1516', '200725_img1516', '200726_img1516', '200727_img1516', '200728_img1516', '200729_img1516', '200730_img1516', '200731_img1516', '200801_img1516', '200802_img1516', '200803_img1516', '200804_img1516', '200805_img1516', '200806_img1516', '200807_img1516', '200808_img1516'};
%days = {'201121_img1520', '201122_img1520', '201123_img1520', '201124_img1520', '201125_img1520', '201126_img1520', '201127_img1520', '201128_img1520', '201129_img1520', '201130_img1520', '201201_img1520', '201202_img1520', '201203_img1520', '201204_img1520', '201205_img1520', '201206_img1520', '201207_img1520', '201208_img1520', '201209_img1520'};

days = {'191118_img1507','200124_img1511','200205_img1510','200331_img1512','200409_img1513','200808_img1516','201209_img1520'};
