function loadPhyOutput
SS_type = input('1 file or 2?');
disp('sample rate assumed to be 30kHz');
disp('pick KS output directory');

folderPath = uigetdir;
% get spikes structure 
temp_spikes = loadKSdir(folderPath);

% specify which ones are MUA 
spikes.MUA = temp_spikes.cgs(temp_spikes.cgs<3)==1;
disp([num2str(numel(spikes.MUA)) ' units; of which ',num2str(sum(spikes.MUA)),' are MUA']);

% keep ones that are not unsorted
keep_clu_id = temp_spikes.cids(temp_spikes.cgs<3);
keep_clu_idx = keep_clu_id + 1;

% find each template ID and keep only ones that are sorted
tempPerClu = findTempForEachClu(temp_spikes.clu, temp_spikes.spikeTemplates);
keep_templates_id = tempPerClu(keep_clu_idx);
keep_templates_idx = keep_templates_id + 1;

% find spike times and ISI violations
spikes.spiketime_S = arrayfun(@(x) double(temp_spikes.st(temp_spikes.clu==x)),keep_clu_id,'un',0)';
spikes.isi_ms = cellfun(@(x) diff(x), spikes.spiketime_S,'un',0);
spikes.isi_violation = cell2mat(cellfun(@(x) sum(x<=0.001)/numel(x),spikes.isi_ms,'un',0));

% find template waveforms and depths 
[~,temp_spikes.spikeDepths,temp_spikes.templateDepths,temp_spikes.templateAmps,~,temp_spikes.templateDuration,temp_spikes.waveform] = templatePositionsAmplitudes(temp_spikes.temps,temp_spikes.winv,temp_spikes.ycoords,temp_spikes.spikeTemplates,temp_spikes.tempScalingAmps);

% keep only ones that are sorted
spikes.template_big = temp_spikes.waveform(keep_templates_idx,:)';
spikes.template_pos = temp_spikes.templateDepths(keep_templates_idx);
spikes.wave_stats.pt_ms = temp_spikes.templateDuration(keep_templates_idx);
spikes.wave_stats.pt_distance = temp_spikes.templateAmps(keep_templates_idx);

tempSplit = strsplit(folderPath,'\');
% tempSplit = strsplit(folderPath,'/');
if SS_type==1
    dFN='1';
else
    dFN='1-2';
    end
saveName = fullfile(folderPath,'..','..');
tempSplit2=strsplit(saveName,'\');
savename=[tempSplit2{1},'\',tempSplit2{2},'\',tempSplit2{3},'\',tempSplit{4},'\',tempSplit2{5},'\',tempSplit2{6},'\'];
filename=[tempSplit{6},'_','datafile00',dFN,'_spikes.mat'];
save([savename,filename],'spikes');
end

