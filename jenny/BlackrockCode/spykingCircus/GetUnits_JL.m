function [spikes] = GetUnits_JL
%% Try to read the waveform information of each template

folderPath = uigetdir;
ver = input('Which version?');

folderIdx = strfind(folderPath,'\');
folderName = folderPath(folderIdx(end)+1:end);
tmpfile = [folderPath,'\',folderName,'.templates-',num2str(ver),'.hdf5'];
probtype = input('which probe type? 1 = miao_poly; 2 = H4_linear \n');

if probtype == 1
probfile = 'Z:\home\jen\Misc\spykingCircus\miao_poly.prb';
else
probfile = 'Z:\home\jen\Misc\spykingCircus\H4_linear.prb';
end
%probfile = 'R:\home\robin\ephys_analysis\spykingCircus\miao_shank.prb';
resultfile = [folderPath,'\',folderName,'.result-',num2str(ver),'.hdf5'];
clusterfile = [folderPath,'\',folderName,'.clusters-',num2str(ver),'.hdf5'];
Fs = 30000;

% in the template file, the first half is the average waveform, the second
% half is the direction of largest variance that is orthogonal to this
% average waveform;

templates_size = double(h5read(tmpfile, '/temp_shape'));
N_e = templates_size(1); % number of channels
N_t = templates_size(2); % the temporal width of the template
temp_x = double(h5read(tmpfile, '/temp_x') + 1);
temp_y = double(h5read(tmpfile, '/temp_y') + 1);
temp_z = double(h5read(tmpfile, '/temp_data'));
templates = sparse(temp_x, temp_y, temp_z, templates_size(1)*templates_size(2), templates_size(3));
templates_size = [templates_size(1) templates_size(2) templates_size(3)/2];
N_temp = templates_size(3); % number of template that is classified, the other half is just the
% get the prob file
pdata = importdata(probfile);
ch_line = strsplit(cell2mat(pdata(cellfun(@(IDX) ~isempty(IDX), strfind(pdata, 'total_nb_channels')))),'= ');
nchannels = str2num(ch_line{2}); % number of total channels
% Extracting channel geometry
ch_info = nan(nchannels,3);
chinfo_lines = find(cellfun(@(IDX) ~isempty(IDX), strfind(pdata, ': ['))); % the first line is not channle infor,
%%
for ch = 1:nchannels
    temp = strsplit(pdata{chinfo_lines(ch+1)}, {' ',':','[',',',']'});
    chIDX = str2num(temp{2})+1;
    ch_info(chIDX,1) = str2num(temp{2})+1; % channel ID
    ch_info(chIDX,2) = str2num(temp{3}); % channel x position [0] superficial, [750] deep
    ch_info(chIDX,3) = str2num(temp{4}); % channel y position
end

channelposlist = ch_info(:,2:3);


Temp = [];
Temp_pos = [];
Temp_big = []; % get one representative template for each cluster
chID = [];
for i_temp = 1:N_temp
    
    Temp(:,:,i_temp) = full(reshape(templates(:, i_temp),[N_t, N_e])); % store the template for each unit
    % find the position of each template
    [Temp_pos(i_temp,1:2),Temp_big(:,i_temp),chID(i_temp,:)]= WFpositioner(Temp(:,:,i_temp),channelposlist);
    %get the spike statistics
    [pt_ms(i_temp,1),pt_ratio(i_temp,1),pt_distance(i_temp,1),asym(i_temp,1),hfw(i_temp,1),endslope(i_temp,1)]=waveStats(Temp_big(:,i_temp));
    
end

%  figure
% plot(Temp_big(:,1))
% hold on
% plot(Temp_big(:,4))

%% Get the spiketimes from a cell and compute the ISI. Spike times is in sampling rate...



for  i_temp = 1:N_temp
    
    tempStr = [];
    
    tempStr = ['/spiketimes/temp_' num2str(i_temp-1)]; % in the hdf5 file, the temp_0 is the first template
    spikes.spiketime_S {i_temp,1}  = double(h5read(resultfile, tempStr))./Fs; % spiketimes in seconds store in the same order as the template
    spikes.isi_ms{i_temp,1} = diff(spikes.spiketime_S {i_temp,1})*1000;
    spikes.isi_vilation(i_temp) = sum(spikes.isi_ms{i_temp,1}<=2)./length( spikes.isi_ms{i_temp,1});% calculate isi violations
    % figure
    % histogram(spikes.isi_ms{i_temp,1},[0:0.3:35],'FaceColor',[0 0 0],'EdgeColor',[0 0 0])
end
%% get the cluster information

for  i_temp = 1:N_temp
    clusterStr = [];
    clusterStr = ['/clusters/temp_' num2str(i_temp-1)];
    spikes.cluster{i_temp,1} = double(h5read(clusterfile, clusterStr)); % clusters store the clusters for multiple dimensions, meaningfull if two clusters are identified within one site
   
%     scatter (spikes.cluster{i_temp,1}(2,:), spikes.cluster{i_temp,1}(3,:))
%     hold on
end
%% store the data and get wave statistics
spikes.channelposlist = channelposlist; % channel position
spikes.template = Temp; % raw templates
spikes.template_pos = Temp_pos; % template positions
spikes.template_big = Temp_big; % each template is stored in each roluw
spikes.chID = chID; % the channel information for each template
spikes.waveform.unit = 'uV';
spikes.wave_stats.pt_ms = pt_ms;
spikes.wave_stats.pt_ratio = pt_ratio;
spikes.wave_stats.pt_distance = pt_distance;
spikes.wave_stats.asym =asym;
spikes.wave_stats.hfw =hfw;
spikes.wave_stats.endslope = endslope;

spikes.threshold = folderName;
%% save file in same folder

save([folderPath(1:folderIdx(numel(folderIdx))) ...
    folderPath(folderIdx(6)+1:folderIdx(7)-1) '_' folderPath(folderIdx(7)+1:folderIdx(8)-1) '_spikes'],'spikes');
end
