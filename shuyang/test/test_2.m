% cross-correlogram:
% use CS as references, SS are targets. for each CS: look at time bin of
% ...ms around it, how many SS in each time bin? average across CS and plot
% a bar graph

analysis_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
channel = '191216_Ch6_MM';
PC = load([analysis_dest channel '.mat']);
allspk_times = PC.spk_time; % in seconds

% channel_CS = '191206_Chan16_filtered';
% spk2_out_all = converted_ns6_191206_SJ_Ch401;
% values_all = spk2_out_all.values;

% times_CS_manul = [0.2297, 0.4851, 0.6406, 1.3295, 1.5270, 1.7875, 1.9705,...
%     2.4067, 2.5029, 2.5157, 2.5316, 2.5442, 2.7046, 2.8440, 3.5368, 3.7003,...
%     4.0935, 4.1293, 4.3670, 4.7245, 5.0406, 5.0882, 6.2229, 6.4445, 7.5780,...
%     7.6047, 7.6876, 7.7412, 9.3366, 10.5376, 10.7435, 10.9578, 13.9465, ...
%     14.2684, 14.3227, 19.7984, 49.1163];% i manully picked these times in spike2
times_CS_manul = [3.5981, 3.7003, 4.0935, 4.3670, 4.7245, 4.7861, 4.8368,...
    7.6046, 7.7758, 8.1574, 10.5932, 12.9371, 14.2684, 14.4404, 21.6858,...
    23.3221, 24.3859, 24.8617, 26.7683, 27.0537, 27.6244, 28.8179, ...
    30.3940, 30.7174, 31.1780, 31.2953, 31.7188, 33.5743, 34.6404,...
    35.6644, 36.6086, 36.7418, 38.7382, 45.8482, 46.6722, 48.0324]; % in seconds

% times_CS_manul = [23.3221, 24.3859, 24.8617, 26.7683, 27.0537, 27.6244, ...
%     28.8179, 30.3940, 30.7174, 31.1780, 31.2953, 31.7188, 33.5743, 34.6404,...
%     35.6644, 36.6086, 36.7418, 38.7382, 45.8482, 46.6722, 48.0324];

% look at 100ms? 
taround = 100/1000; % in seconds
corr_mat = zeros(length(times_CS_manul),20);


for i = 1: length(times_CS_manul)
    t1 = times_CS_manul(i) - taround;
    %t2 = times_CS_manul(i) + taround;
    %bin in every 10ms
    for t = 1:20
        corr_mat(i,t) = length(find(allspk_times > t1+(t-1)*10/1000 & allspk_times <= t1+t*10/1000));
        a = t1+(t-1)*10/1000
        b = t1+t*10/1000
    end
end

ave_corr_mat = mean(corr_mat,1);
crosscorrelogram = figure;
x = (0:0.01:0.2);
histogram('BinEdges',x, 'BinCounts', ave_corr_mat);
vline(0.1,'r');
title('cross correlogram');
ylabel('number of simple spikes');xlabel('time(s)');



