function param = spikeWaveform_separator(spikes)

in = 'n';
while ~isempty(in)
figure
plot(spikes.wave_stats.pt_distance,spikes.wave_stats.pt_ms,'ro');
fix_axes(gcf,15,'peak-trough distance','peak-trough time');
[~,y] = ginput;
hold on
IN = spikes.wave_stats.pt_ms > y;
plot(spikes.wave_stats.pt_distance(IN) ,spikes.wave_stats.pt_ms(IN), 'ko');

figure
hold on
timeVector = make_time(spikes.template_big(:,IN),30000,1);
plot(timeVector,spikes.template_big(:,IN),'k');
plot(timeVector,spikes.template_big(:,~IN),'r');


in = input('proceed?');
end
close all

param = IN;
end