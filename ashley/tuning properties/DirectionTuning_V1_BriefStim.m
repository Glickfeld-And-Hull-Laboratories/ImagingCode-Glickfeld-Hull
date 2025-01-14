SubNum = '613';
mouse = 'AW13';
date = '150511';
time = '1523';
ImgFolder = '006';

% analysis folder
try
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(date,ImgFolder)
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
end

%% Parameters
nON = (input.nScansOn);
nOFF = (input.nScansOff);
nStim = input.gratingDirectionStepN;

%% downsample and register data

%remove negative data (by addition)
data_sub = data-min(min(min(data,[],1),[],2),[],3);
clear data

% register
DirFolder = '005';
fileDirMasks = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, DirFolder);
cd(fileDirMasks);
load('regImg.mat')


% data_avg = mean(data_sub(:,:,101:200),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%% neuropil subtract
data_TC = stackGetTimeCourses(data_reg,mask_cell);
figure; tcOffsetPlot(data_TC)
nCells = size(data_TC,2);
buf = 4;
np = 6;
neuropil = imCellNeuropil(mask_cell,buf,np);
% for i = 1:nCells
%     npTC(:,i) = stackGetTimeCourses(data_reg,neuropil(:,:,i));
% end
npTC = stackGetTimeCourses(data_reg,neuropil);

%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_TC-tcRemoveDC(npTC*ii(i)));
end
[max_skew ind] =  max(x,[],1);
skew(buf,:) = max_skew;
np_w = 0.01*ind;
npSubTC = data_TC-bsxfun(@times,tcRemoveDC(npTC),np_w);

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
save('neuropil.mat','neuropil','npTC','npSubTC');

data_TC = npSubTC;
%%
% nRep = size(data_reg,3)./((nON+nOFF)*nStim);
nRep = size(data_TC,1)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
%% create dF/F stack
data_reg = double(data_reg);

nOFF_ind = zeros(1,(nOFF*nStim*nRep));
start = 1;
for iStim = 1:(nRep*nStim)
    nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
    start = start+nOFF;
end

nON_ind = setdiff(1:size(data_reg,3),nOFF_ind);
nON_avg = mean(data_reg(:,:,nON_ind),3);
nOFF_avg = mean(data_reg(:,:,nOFF_ind),3);

% dF average F
dF_data = bsxfun(@minus,data_reg, nOFF_avg);
dFoverF = bsxfun(@rdivide,dF_data,nOFF_avg);

max_dF = max(dF_data,[],3);
maxDFoverF = max(dFoverF,[],3);
figure; imagesq(max_dF); colormap(gray)

% dirTuningFolder = '005';
% fileDir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirTuningFolder);
% cd(fileDir);
load('mask&TCDir.mat')

% bwout = imCellEditInteractive(max_dF);
% mask_cell = bwlabel(bwout);

data_TC = stackGetTimeCourses(data_reg,mask_cell);
figure; tcOffsetPlot(data_TC)

%%
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
save('dataTC.mat','data_TC');

%%
nTrials = (nStim.*nRep);
DirectionDeg = cell2mat(input.tGratingDirectionDeg);
Dirs = unique(DirectionDeg);

%% dF/F by trial
% F per trial
stimOFF_ind = 1:nOFF+nON:size(data_TC,1);

dF_data = zeros(size(data_TC));
dFoverF_data = zeros(size(data_TC));
for i = 1:nTrials
    indAll = stimOFF_ind(i):stimOFF_ind(i)+(nOFF+nON-1);
    indF = stimOFF_ind(i)+5:stimOFF_ind(i)+(nOFF-1);
    dF_data(indAll,:) = bsxfun(@minus,data_TC(indAll,:),mean(data_TC(indF,:),1));
    dFoverF_data(indAll,:) = bsxfun(@rdivide,dF_data(indAll,:),mean(data_TC(indF,:),1));
end
%% dF/F (by cell) for each stimulus type

% find on indices for the first frame of each stimulus start period and iti (Off) period

stimON_ind = nOFF+1:nOFF+nON:size(dFoverF_data,1);

% dFoverFCellsTrials = zeros(nON+nOFF,size(dFoverF_data,2),nTrials);
% for i = 1:nTrials
%     dFoverFCellsTrials(:,:,i) = dFoverF_data(stimOFF_ind(i):stimOFF_ind(i)+nON+nOFF-1,:);
% end
% 
% dFoverF_meanDirResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim);
% for i = 1:nStim
%     trials = find(DirectionDeg == Dirs(i));
%     dFoverF_meanDirResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
% end

% sort data_TC into 20 frame (10 pre, 10 post) traces around stimON 

dFoverFCellsTrials = zeros(20+nON,size(dFoverF_data,2),nTrials-1);
for i = 1:nTrials-1
    dFoverFCellsTrials(:,:,i) = dFoverF_data(stimON_ind(i)-10:stimON_ind(i)+(nON-1)+10,:);
end

dFoverF_meanDirResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim);
for i = 1:nStim
    trials = find(DirectionDeg(:,1:end-1) == Dirs(i));
    dFoverF_meanDirResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
end

figure;
for i = 1:nStim
    plot(dFoverF_meanDirResp(:,158,i));
    hold on
    vline(10,'k');
    hold on
end

%% find magnitude of response to stim

dFoverF_meanOFFDirResp = (squeeze(mean(dFoverF_meanDirResp(1:10,:,:),1)));

dFoverF_meanONDirResp = (squeeze(mean(dFoverF_meanDirResp(11:end,:,:),1)));
DirRespPerCell = dFoverF_meanONDirResp;
% DirRespPerCell = dFoverF_meanONDirResp./dFoverF_meanOFFDirResp;

figure;
plot(DirRespPerCell(158,:))

%% find direction preference

DirRespPerCell_sub = DirRespPerCell-min(min(min(DirRespPerCell,[],1),[],2));

[dirPref_val,dirPref_ind] = max(DirRespPerCell_sub,[],2);

pref2orth_dir = [nStim/2+1:nStim 1:nStim/2];
dirOrth_ind = zeros(size(dirPref_ind));
for i = 1:nStim
    cells = find(dirPref_ind == i);
    dirOrth_ind(cells) = pref2orth_dir(i);
end

dirOrth_val = zeros(size(dirPref_ind)); 
for i = 1:size(data_TC,2)
    ind = dirOrth_ind(i);
    cell = DirRespPerCell_sub(i,:);
    dirOrth_val(i) = cell(ind);
end

cellDSI = zeros(size(dirPref_ind));
cellDSI = (dirPref_val - dirOrth_val)./(dirPref_val + dirOrth_val);
dirSlctvCells = find(cellDSI > 0.3);

%% find orientation preference
dFoverF_meanOriResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim/2);
for i = 1:nStim/2
    trials = find(DirectionDeg(1:end-1,:) == Dirs(i) | DirectionDeg(1:end-1,:) == Dirs(i+nStim/2));
    dFoverF_meanOriResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
end

dFoverF_meanOFFOriResp = (squeeze(mean(dFoverF_meanOriResp(1:10,:,:),1)))+1;

dFoverF_meanONOriResp = (squeeze(mean(dFoverF_meanOriResp(11:end,:,:),1)))+1;

OriRespPerCell = dFoverF_meanONOriResp;

figure;
plot(OriRespPerCell(10,:))

OriRespPerCell_sub = OriRespPerCell-min(min(min(OriRespPerCell,[],1),[],2));
[oriPref_val,oriPref_ind] = max(OriRespPerCell_sub,[],2);

pref2orth_ori = [(nStim/2)/2+1:nStim/2 1:(nStim/2)/2];
oriOrth_ind = zeros(size(oriPref_ind));
for i = 1:nStim/2
    cells = find(oriPref_ind == i);
    oriOrth_ind(cells) = pref2orth_ori(i);
end

oriOrth_val = zeros(size(oriPref_ind)); 
for i = 1:size(data_TC,2)
    ind = oriOrth_ind(i);
    cell = OriRespPerCell_sub(i,:);
    oriOrth_val(i) = cell(ind);
end

cellOSI = zeros(size(oriPref_ind));
cellOSI = (oriPref_val - oriOrth_val)./(oriPref_val + oriOrth_val);
oriSlctvCells = find(cellOSI > 0.3);

save('TuningPreferences.mat','oriPref_ind','dirPref_ind','dirSlctvCells','oriSlctvCells')

%% 
cMap = colormap(jet(nStim));
start = 1;
errbar = zeros(nStim,size(data_TC,2));
for ifig = 1:ceil(size(data_TC,2)/15)
figure;
for iplot = 1:15
%     cell = baselineStimRespIndex_V(start);
    cell = start;
    ymax = .1;
    subplot(4,4,iplot);
    for i = 1:nStim
        plot(dFoverF_meanDirResp(:,cell,i),'color',cMap(i,:));
        hold on
        vline(10,'k');
        hold on
        title(['Cell ' num2str(cell)]);
        hold on
        ymax_i = max(dFoverF_meanDirResp(:,cell,i),[],1);
        if ymax_i > ymax
            ymax = ymax_i;
        end
        axis([0 30 -0.05 ymax]);
        hold on
        legendInfo{i} = [num2str(Dirs(i)) ' degrees'];
        hold on

    end   
    if iplot == 1
        legend(legendInfo,'Location','SouthEast')
    end
    start = start+1;
    hold on
%     if ismember(cell,baselineStimRespIndex_V) > 0
%         set(subplot(4,4,iplot),'color',[0.9 0.9 0.9])
%     end
end
end

dFoverFDirResp = zeros(nStim,size(data_TC,2));
errbar = zeros(nStim,size(data_TC,2));
for i = 1:nStim 
    trials = find(DirectionDeg(:,1:end-1) == Dirs(i));
    dFoverFDirResp(i,:) = squeeze(mean(mean(dFoverFCellsTrials(end-10:end,:,trials),1),3));
    errbar(i,:) = std(mean(dFoverFCellsTrials(11:16,:,trials),1),[],3)/sqrt(size(dFoverFCellsTrials(11:16,:,trials),3));
end


start = 1;
for ifig = 1:ceil(size(data_TC,2)/15)
figure;
for iplot = 1:15
%     cell = baselineStimRespIndex_V(start);
    cell = start;
    ymax = .1;
    subplot(4,4,iplot);
    errorbar(dFoverFDirResp(:,cell),errbar(:,cell),'k');
    hold on
    ymax_i = max(dFoverFDirResp(:,cell),[],1);
    if ymax_i > ymax
        ymax = ymax_i;
    end
    title(['Cell ' num2str(cell)]);
    hold on
    axis([0 length(Dirs) -0.05 ymax]);
    hold on
%     if ismember(cell,baselineStimRespIndex_V) > 0
%         set(subplot(4,4,iplot),'color',[0.9 0.9 0.9])
%     end
    hold on
    start = start+1;
end
end

