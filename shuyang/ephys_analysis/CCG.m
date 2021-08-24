% cross-correlogram:
% use CS as references, SS are targets. for each CS: look at time bin of
% ...ms around it, how many SS in each time bin? average across CS and plot
% a bar graph
clear;
analysis_dest = 'Y:\home\nathan\DATA\S-probe\Awake\useful_ephys\analysis\';
channel = '191213_D002_Ch17';
PC = load([analysis_dest channel '.mat']);
allspk_times = PC.spk_time; % in seconds
CS_negPinx = PC.CS_negPinx; % in 30000Hz
CS_negPinx = CS_negPinx/30000; %in seconds

% look at 100ms? 
taround = 100/1000; % in seconds
corr_mat = zeros(length(CS_negPinx),20);

for i = 1: length(CS_negPinx)
    t1 = CS_negPinx(i) - taround;
    %bin in every 5ms
    for t = 1:40
        corr_mat(i,t) = length(find(allspk_times > t1+(t-1)*5/1000 & allspk_times <= t1+t*5/1000));
        a = t1+(t-1)*10/1000
        b = t1+t*10/1000
    end
end

ave_corr_mat = mean(corr_mat,1);
crosscorrelogram = figure;
x = (0:0.005:0.2);
histogram('BinEdges',x, 'BinCounts', ave_corr_mat);
vline(0.1,'r');
title(['cross correlogram' channel]);
ylabel('number of simple spikes');xlabel('time(s)');
savefig([analysis_dest 'figures\' channel '_CCG']);



