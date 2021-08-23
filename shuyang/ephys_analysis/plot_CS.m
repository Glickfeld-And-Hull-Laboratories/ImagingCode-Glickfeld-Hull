%look at all spikes
%load data, drag the file into matlab
channel = '191213_D002_Ch17';
% drag filtered data in matlab
spk2_out_all = converted_ns6_191213_D002_SJ_Ch403;
analysis_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
values_all = spk2_out_all.values;

times_CS_manul = [12.9401, 14.9819, 15.1989, 15.5197, 16.3977, 18.1462,...
    18.5976, 21.9666, 22.5485, 22.9457];% this is what Court manually picked 
values_inx = times_CS_manul*30000; % we're recroding at 30kHz
values_inx_mat = zeros(0.01*30000,length(values_inx));% 10 ms = 0.01s, 0.01*30000 is how long you want
for i = 1:length(times_CS_manul)
    values_inx_mat(:,i) = values_inx(i):values_inx(i)+0.01*30000-1;
end
values_inx_mat = int32(values_inx_mat);% convert double to integer for indexing
values_CS_mat = values_all(values_inx_mat);
% find the first negative peak and align these events to the negative peak
x1 = 1:1:length(values_CS_mat);
CS_negPinx = zeros(1,length(times_CS_manul));
for i = 1:length(times_CS_manul)
    [p,k] = findpeaks(-(values_CS_mat(:,i)));
    %figure; plot(x1,values_CS_mat(:,i)); hold on; plot(k(1),-(p(1)),'ro');hold off;% test if the first peak is the real negative peak, if not, adjust time manul time points above
    CS_negPinx(i) = values_inx_mat(k(1),i);
end

CS_inx_aligned_mat = zeros(0.01*30000,length(times_CS_manul));% 10 ms = 0.01s, 0.01*30000 is how long you want
for i = 1:length(times_CS_manul)
    CS_inx_aligned_mat(:,i) = CS_negPinx(i)-0.001*30000:CS_negPinx(i)+0.009*30000-1;
end
CS_inx_aligned_mat = int32(CS_inx_aligned_mat);
values_CS_aligned_mat = values_all(CS_inx_aligned_mat);
ave_CS_aligned = mean(values_CS_aligned_mat,2);
save([analysis_dest channel '.mat'],'values_CS_aligned_mat','ave_CS_aligned',...
    'times_CS_manul','CS_inx_aligned_mat','CS_negPinx','-append');


CS_picked = figure;
x = (1:0.01*30000)*1000/30000;
xmat = repmat(x,size(values_CS_aligned_mat,2),1); xmat = xmat';
plot(xmat,values_CS_aligned_mat,'color',[0.8 0.8 0.8]); 
ylabel('Volts');
xlabel('time (ms)');
savefig([analysis_dest 'figures\191213_D002_Ch17_CS_indi']);
figure; plot(x,ave_CS_aligned,'color','k'); 
ylabel('Volts');
xlabel('time (ms)');
savefig([analysis_dest 'figures\191213_D002_Ch17_CS_ave']);



% fig_name = [channel '_allspk'];
% path = [analysis_dest 'figures\'];
% %savefig(allspike, [path,fig_name]);

% times_CS_manul = [0.2297, 0.4851, 0.6406, 1.3295, 1.5270, 1.7875, 1.9705,...
%     2.4067, 2.5029, 2.5157, 2.5316, 2.5442, 2.7046, 2.8440, 3.5368, 3.7003,...
%     4.0935, 4.1293, 4.3670, 4.7245, 5.0406, 5.0882, 6.2229, 6.4445, 7.5780,...
%     7.6047, 7.6876, 7.7412, 9.3366, 10.5376, 10.7435, 10.9578, 13.9465, ...
%     14.2684, 14.3227, 19.7984, 49.1163];% i manully picked these times in spike2
% times_CS_manul = [3.5981, 3.7003, 4.0935, 4.3670, 4.7245, 4.7861, 4.8368,...
%     7.6046, 7.7758, 8.1574, 10.5932, 12.9371, 14.2684, 14.4404, 21.6858,...
%     23.3221, 24.3859, 24.8617, 26.7683, 27.0537, 27.6244, 28.8179, ...
%     30.3940, 30.7174, 31.1780, 31.2953, 31.7188, 33.5743, 34.6404,...
%     35.6644, 36.6086, 36.7418, 38.7382, 45.8482, 46.6722, 48.0324];
% times_CS_manul = [23.3221, 24.3859, 24.8617, 26.7683, 27.0537, 27.6244, ...
%     28.8179, 30.3940, 30.7174, 31.1780, 31.2953, 31.7188, 33.5743, 34.6404,...
%     35.6644, 36.6086, 36.7418, 38.7382, 45.8482, 46.6722, 48.0324];
