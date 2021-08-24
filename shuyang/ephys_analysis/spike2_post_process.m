%look at all spikes
%load data, drag the file into matlab

channel = '191206_Ch16';
spk2_out_all = converted_ns6_191206_SJ_Ch37;
codes_all = spk2_out_all.codes;
codes_all = codes_all(:,1);
values_all = spk2_out_all.values;
times_all = spk2_out_all.times;
ntotal = spk2_out_all.length;
vectotal = 1:ntotal;
spks_inx = vectotal(codes_all==1); %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pay attention to the codes
spks_inx = sort(spks_inx);% make all indexes in ascending order, otherwise this will be a problem if codes_all == more than 1 number
spk_time = times_all(spks_inx);
all_spk = values_all(spks_inx,:);
ISI = diff(spk_time)*1000;%second to ms

analysis_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
save([analysis_dest channel '.mat'],'codes_all','values_all','times_all',...
    'spks_inx','spk_time','all_spk','ISI');

%plot ISI histogram
ISIH = figure; 
histogram(ISI(ISI<1000),'BinWidth',5);
xlabel('time(ms)');
ylabel('counts/bin');
title('ISIH-191206-Ch16'); %!!!!!!!!!!!!!pay attention to the channel number
fig_name = '191206_Ch16_ISIH';%!!!!!!!!!!!!!pay attention to the channel number
path = [analysis_dest 'figures\'];
savefig(ISIH, [path,fig_name]);

%plot all spikes
allspike = figure;
plot(all_spk','color',[0.8 0.8 0.8]); hold on;
plot(mean(all_spk),'color','k');
ylabel('Volts');
title('all spikes,191206-Ch16'); %!!!!!!!!!!!!!pay attention to the channel number
fig_name = [channel '_allspk'];
path = [analysis_dest 'figures\'];
%savefig(allspike, [path,fig_name]);

%{
%load data, drag the file into matlab
spk2_out_CS = converted_ns6_191213_D002_SJ_Ch39;
codes_CS = spk2_out_CS.codes;
values_CS = spk2_out_CS.values;
nCStotal = spk2_out_CS.length;
times_CS = spk2_out_CS.times;
vecCStotal = 1:nCStotal;
CS_inx = vecCStotal(codes_CS==1);
CS_time = times_CS(codes_CS==1);
cspk = values_CS(CS_inx,:);
% complexspike = figure;
% plot(cspk','color',[0.8 0.8 0.8]);hold on;
% plot(mean(cspk),'color','k');
% ylim([-0.3 0.2]);
% xlabel('time(ms)');
% ylabel('Volts');
% title('Cspk, 191213-D002-Ch16');
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontSize',18);
% fig_name = '191213_D002_Ch16_CS';
% path = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
% savefig(complexspike, [path,fig_name]);


spk2_out_all = converted_ns6_191213_D002_SJ_Ch38;
codes_all = spk2_out_all.codes;
values_all = spk2_out_all.values;
times_all = spk2_out_all.times;
ntotal = spk2_out_all.length;
vectotal = 1:ntotal;
all_spk_inx = vectotal(codes_all==1);
%all_spikes = values_all(all_spk_inx,:);
all_spikes_time = times_all(all_spk_inx);
SS_time = setdiff(all_spikes_time,CS_time);
SS_inx = ismember(all_spikes_time,SS_time);
SS = values_all(SS_inx,:);
simplespike = figure;
plot(SS','color',[0.8 0.8 0.8]);hold on;
plot(mean(SS),'color','k');
ylim([-0.4 0.3]);
xlabel('time(ms)');
ylabel('Volts');
title('SS, 191213-D002-Ch42');
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',18);
fig_name = '191213_D002_Ch42_SS';
path = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
savefig(simplespike, [path,fig_name]);
%}
