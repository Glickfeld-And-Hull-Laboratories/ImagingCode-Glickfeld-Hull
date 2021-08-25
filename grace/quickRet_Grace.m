% for doing a retinotopy when using multiple positions on one screen
clear all
clear global
%% 
date = '210424';
mouse = 'i1349';
ImgFolder = '001';
time = '1529';
run = '000';
doReg = 0;
nrun = size(ImgFolder,1);

% tHostname = lower(hostname);
% [s,tUsername] = dos('ECHO %USERNAME%');

run_str = ['runs-' ImgFolder(1,:)];
if nrun>1
    run_str = [run_str '-' ImgFolder(nrun,:)];
end

% if strcmp(tUsername(1:5),'lindsey')
CD = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\data\2P_images\' mouse '\' date '\' ImgFolder];
% else
% error('Not Grace')    
% end
cd(CD);
imgMatFile = [ImgFolder '_000_' run '.mat'];
load(imgMatFile);

nframes = info.config.frames;
data = sbxread([ImgFolder '_000_' run],0,nframes);
fprintf(['Loaded ' num2str(nframes) ' frames \r\n'])

%default is to use only green PMT if there are two
if size(data,1) > 1
    data = data(1,:,:,:);
end
data = squeeze(data);
   
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
load(fName);

nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = size(input.tGratingDirectionDeg,2);
for irun = 1:nrun
    temp(irun) = input;
    if isfield(input, 'nScansOn')
        nOn = temp(irun).nScansOn;
        nOff = temp(irun).nScansOff;
        ntrials = size(temp(irun).tGratingDirectionDeg,2);
     if nframes<ntrials*(nOn+nOff)
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
     end
    end
end
input = concatenateDataBlocks(temp);
ntrials = size(input.tGratingDirectionDeg,2);

t = 1000;
nep = floor(size(data,3)./t);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*t):200+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); end

%% 
if doReg
data_avg = mean(data(:,:,1000:1500),3);
[out data_reg] = stackRegister(data,data_avg);
data = data_reg;
clear data_reg
end
    
sz = size(data);
data = data(:,:,1:(nOn+nOff)*ntrials);
data_tr = reshape(data,[sz(1), sz(2), nOn+nOff, ntrials]);
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
data_df = bsxfun(@minus, double(permute(data_tr,[1 2 4 3])), squeeze(data_f)); 
data_dfof = permute(bsxfun(@rdivide,data_df, squeeze(data_f)),[1 2 4 3]); 
clear data_f data_df data_tr

Az = celleqel2mat_padded(input.tGratingAzimuthDeg);
El = celleqel2mat_padded(input.tGratingElevationDeg);
Azs = unique(Az);
Els = unique(El);
if min(Els,[],2)<0
    Els = fliplr(Els);
end
nStim = length(Azs).*length(Els);
Stims = [];
data_dfof_avg = zeros(sz(1), sz(2), nOn+nOff, length(Azs).*length(Els));
start = 1;
for iEl = 1:length(Els)
    ind1 = find(El == Els(iEl));
    for iAz = 1:length(Azs)
        Stims = [Stims; Els(iEl) Azs(iAz)];
        ind2 = find(Az == Azs(iAz));
        ind = intersect(ind1,ind2);
        data_dfof_avg(:,:,:,start) = mean(data_dfof(:,:,:,ind),4);
        start = start +1;
    end
end
clear data_dfof
myfilter = fspecial('gaussian',[20 20], 0.5);
data_dfof_avg_all = squeeze(mean(imfilter(data_dfof_avg(:,:,nOff:nOff+nOn,:),myfilter),3));
data_dfof_max = max(data_dfof_avg_all,[],3);
figure;imagesc(data_dfof_max)
% if strcmp(tUsername(1:5),'grace')
    mkdir(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder]);
    print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder '\' date '_' mouse '_' ImgFolder '_data_dfof_max.pdf'], '-dpdf','-bestfit')
% end

img_avg_resp = zeros(1,nStim);
figure; 
for i = 1:nStim
    subplot(length(Els),length(Azs),i); 
    imagesc(data_dfof_avg_all(:,:,i)); 
    colormap(gray)
    title(num2str(Stims(i,:)))
    img_avg_resp(i) = mean(mean(mean(data_dfof_avg_all(:,:,i),3),2),1);
    %clim([0 max(data_dfof_avg_all(:))./2])
end
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder '\' date '_' mouse '_' ImgFolder '_retinotopy.pdf'], '-dpdf','-bestfit')


figure;
heatmap = imagesc(fliplr(rot90(reshape(img_avg_resp,length(Els),length(Azs)),3)));
heatmap.Parent.YTick = 1:length(Els);
heatmap.Parent.YTickLabel = strread(num2str(Els),'%s');    
heatmap.Parent.XTick = 1:length(Azs);
heatmap.Parent.XTickLabel = strread(num2str(Azs),'%s');
xlabel('Azimuth');
ylabel('Elevation');
colorbar
caxis([-0.1 0.1])
% print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder '\' date '_' mouse '_' ImgFolder '_ret_colorbar.pdf'], '-dpdf','-bestfit')

figure; imagesc(mean(data,3));
axis off
title([mouse ' ' date])
if strcmp(tUsername(1:5),'grace')
% print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P\' date '_' mouse '\' date '_' mouse '_' ImgFolder '\' date '_' mouse '_' ImgFolder '_FOV.pdf'], '-dpdf','-bestfit')
end