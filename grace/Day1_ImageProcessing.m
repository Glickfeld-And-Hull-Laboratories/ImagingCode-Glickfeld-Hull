clear all
clear global
%% get path names
date = '210603';
ImgFolder = strvcat('003');
time = strvcat('1425');
mouse = 'i1802';
run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
ref_str = catRunName(run, size(run,1));
gl_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\2P_Imaging';
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\grace\Analysis\2P';
behav_fn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
%% load and register
data = [];
clear temp
trial_n = [];
offset = 0;
for irun = 1:nrun
    CD = fullfile(gl_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);
    fName = fullfile(behav_fn, ['data-' mouse '-' date '-' time(irun,:) '.mat']);
    load(fName);

    nframes = info.config.frames;
    fprintf(['Reading run ' num2str(irun) '- ' num2str(nframes) ' frames \r\n'])
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    
    temp(irun) = input;
    if isfield(input, 'nScansOn')
        nOn = temp(irun).nScansOn;
        nOff = temp(irun).nScansOff;
        ntrials = size(temp(irun).tGratingDirectionDeg,2);

        data_temp = squeeze(data_temp);
        if nframes>ntrials*(nOn+nOff)
            nframes = ntrials*(nOn+nOff);
            data_temp = data_temp(:,:,1:ntrials*(nOn+nOff));
        elseif nframes<ntrials*(nOn+nOff)
            temp(irun) = trialChopper(temp(irun),1:ceil(nframes./(nOn+nOff)));
        end
    end
    
    offset = offset+nframes;

    data_temp = squeeze(data_temp);
    data = cat(3,data,data_temp);
    trial_n = [trial_n nframes];
end
input = concatenateDataBlocks(temp);
clear data_temp
clear temp

 %% Choose register interval
t = 2000;
nep = floor(size(data,3)./t);
[n n2] = subplotn(nep);
figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data(:,:,1+((i-1)*t):500+((i-1)*t)),3)); title([num2str(1+((i-1)*t)) '-' num2str(500+((i-1)*t))]); end

%% Register data

data_avg = mean(data(:,:,46001:46500),3); 

if exist(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    load(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
    [outs, data_reg]=stackRegister_MA(data(:,:,:),[],[],out);
    clear out outs
else
    [out, data_reg] = stackRegister(data,data_avg);
    data_reg_avg = mean(data_reg,3);
    reg = data_reg_avg;
    mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str]))
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'data_reg_avg', 'out', 'data_avg')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')
end
clear data


%% test stability
% figure; for i = 1:nep; subplot(n,n2,i); imagesc(mean(data_reg(:,:,1+((i-1)*nframes):500+((i-1)*nframes)),3)); title([num2str(1+((i-1)*nframes)) '-' num2str(500+((i-1)*nframes))]); end
figure; imagesq(data_reg_avg); truesize;
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_avg.pdf']),'-dpdf','-bestfit')

%% find activated cells
%max by trial type    
if isfield(input, 'nScansOn')
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    if nOn>29
        sz = size(data_reg);
%         make data_tr smaller for dfof images to not use all of the memory
        data_tr = reshape(single(data_reg(:,:,1:28800)),[sz(1), sz(2), nOn+nOff, 320]);
        data_f = mean(double(data_tr(:,:,nOff/2:nOff,:)),3);
        data_df = bsxfun(@minus, double(data_tr), data_f); 
        data_dfof = bsxfun(@rdivide,data_df, data_f); 
        clear data_tr clear data_f clear data_df
    else
        sz = size(data_reg);
        data_tr = zeros(sz(1),sz(2), 100, ntrials-1);
        for itrial = 1:ntrials-1
            data_tr(:,:,:,itrial) = data_reg(:,:,((itrial-1)*(nOn+nOff))+71:170+((itrial-1)*(nOn+nOff)));
        end
        data_f = mean(data_tr(:,:,1:50,:),3);
        data_df = bsxfun(@minus, double(data_tr), data_f); 
        data_dfof = bsxfun(@rdivide,data_df, data_f); 
        clear data_tr clear data_f clear data_df
    end
end

if input.doDirStim
%     obtaining avg and max dfof images for each dir
    Dir = cell2mat_padded(input.tGratingDirectionDeg);
    Dir = Dir(1:ntrials);
    Dirs = unique(Dir);
    data_dfof_avg = zeros(sz(1),sz(2),length(Dirs));
    nDirs = length(Dirs);
    [n n2] = subplotn(nDirs);
    figure;
    for idir = 1:length(Dirs)
        if nOn>29
            ind = find(Dir(1:160) == Dirs(idir));
        else
            ind = find(Dir(1:ntrials-1) == Dirs(idir));
        end
        data_dfof_avg(:,:,idir) = mean(mean(data_dfof(:,:,nOff+1:nOn+nOff,ind),3),4);
        subplot(n,n2,idir)
        imagesc(data_dfof_avg(:,:,idir))
    end
    clear data_dfof
    myfilter = fspecial('gaussian',[20 20], 0.5);
    data_dfof_avg_all = imfilter(data_dfof_avg,myfilter);
    data_dfof_max = max(data_dfof_avg_all,[],3);
    
%     average dfof image for each direction
    figure; 
    Stims = Dirs;
    nStim = length(Dirs);
    [n n2] = subplotn(nDirs);
    data_dfof_avg_ori = zeros(sz(1), sz(2), nDirs/2);
    for i = 1:nStim 
        subplot(n,n2,i); 
        imagesc(data_dfof_avg_all(:,:,i));
        clim([0 max(data_dfof_avg_all(:))])
        title(num2str(Dirs(i)))
        colormap(gray)
        if i<(nDirs/2)+1
            data_dfof_avg_ori(:,:,i) = mean(data_dfof_avg_all(:,:,[i i+nDirs/2]),3);
        end
    end
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_16Stim.pdf']), '-dpdf')

%     average dfof image for each orientation
    figure;
    [n n2] = subplotn(nDirs/2);
    for i = 1:nStim/2
        subplot(n,n2,i)
        imagesc(data_dfof_avg_ori(:,:,i));
        clim([0 max(data_dfof_avg_ori(:))])
        title(num2str(Dirs(i)))
        axis off
    end
    subplot(n,n2,i+1)
    imagesc(max(data_dfof_avg_ori,[],3))
    title('dfof Max')
    axis off
    data_dfof = cat(3,data_dfof_avg_ori,max(data_dfof_avg_ori,[],3));
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_activeCells.pdf']), '-dpdf')

    figure;
    imagesc(max(data_dfof_avg_ori,[],3))
    title('dfof Max')
    axis off
    print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofMax.pdf']), '-dpdf')
end

 save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimActFOV.mat']), 'data_dfof_max', 'data_dfof_avg_all')


%% cell segmentation 
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1),sz(2));
mask_data = data_dfof_avg_all;

% click on cells in each direction's avg dfof image
for iStim = 1:size(data_dfof_avg_all,3)
  mask_data_temp = mask_data(:,:,iStim);
  mask_data_temp(find(mask_exp >= 1)) = 0;
  bwout = imCellEditInteractive(mask_data_temp);
  mask_all = mask_all+bwout;
  mask_exp = imCellBuffer(mask_all,3)+mask_all;
  close all
end
mask_cell = bwlabel(mask_all);
figure; imagesc(mask_cell)


figure; 
[n n2] = subplotn(nStim);
for i = 1:nStim; 
    subplot(n,n2,i); 
    shade_img = imShade(data_dfof_avg_all(:,:,i), mask_all);
    imagesc(shade_img)
    if input.doSizeStim
    title([num2str(szs(i)) ' deg'])
    elseif input.doRetStim
        title([num2str(Stims(i,:))])
    end
    clim([0 max(data_dfof_avg_all(:))])
    colormap(gray)
end
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOV_overlay.pdf']), '-dpdf')


mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'data_dfof_max', 'mask_cell', 'mask_np')

% creating pixel correlation image from a smaller data reg stack
data_reg_3hz = stackGroupProject(data_reg,5);
pix = getPixelCorrelationImage(data_reg_3hz);
pix(isnan(pix))=0;
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pixel.mat']),'pix')

clear data_dfof data_dfof_avg max_dfof mask_data mask_all mask_2 data_base data_base_dfof data_targ data_targ_dfof data_f data_base2 data_base2_dfof data_dfof_dir_all data_dfof_max data_dfof_targ data_avg data_dfof2_dir data_dfof_dir 
%% neuropil subtraction
data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_tc_down = stackGetTimeCourses(stackGroupProject(data_reg,5), mask_cell);
nCells = size(data_tc,2);
sz = size(data_reg);
down = 5;
data_reg_down  = stackGroupProject(data_reg,down);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf([' Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc(:,:)-bsxfun(@times,tcRemoveDC(np_tc(:,:)),np_w);

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']), 'input')

%% 
ntrials = size(input.tGratingDirectionDeg,2);
Dir = cell2mat_padded(input.tGratingDirectionDeg);
Dir = Dir(1:ntrials);
    Dirs = unique(Dir);
    nDirs = length(Dirs);
if isfield(input, 'nScansOn')
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    nCells = size(npSub_tc,2);
    if nOn>29
        data_mat = zeros(nOn+nOff, nCells, ntrials);
        for itrial = 1:ntrials
            data_mat(:,:,itrial) = npSub_tc(1+((itrial-1).*(nOn+nOff)):(itrial.*(nOn+nOff)),:);
        end
        data_f = mean(data_mat(nOff/2:nOff,:,:),1);
    else
        data_mat = zeros(100, nCells, ntrials-1);
        for itrial = 1:ntrials-1
            data_mat(:,:,itrial) = npSub_tc(((itrial-1)*(nOn+nOff))+71:170+((itrial-1)*(nOn+nOff)),:);
        end
        data_f = mean(data_mat(1:50,:,:),1);
    end
    data_df = bsxfun(@minus, data_mat, data_f);
    data_dfof = bsxfun(@rdivide, data_df, data_f);
    
    ndir = length(Dirs);
    [n, n2] = subplotn(nCells);
    h_dir = zeros(nCells, ndir);
    p_dir = zeros(nCells, ndir);
    base_win = 50:60;
    resp_win = 70:90;
    base = squeeze(mean(data_dfof(base_win,:,:),1));
    resp = squeeze(mean(data_dfof(resp_win,:,:),1));
    dir_resp = zeros(nCells,ndir);
    [x y] = ttest(resp', base', 'tail','right');
    no_match = find(isnan(x));
    max_dir = zeros(nCells,1);
    figure;
    for i = 1:nCells
        if ~sum(no_match==i)
        subplot(n, n2, i)
            for idir = 1:ndir
                if nOn>29
                    ind = find(Dir == Dirs(idir));
                else
                    ind = find(Dir(1:ntrials-1) == Dirs(idir));
                end
                [h_dir(i,idir), p_dir(i,idir)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(ndir-1));
                if h_dir(i,idir)
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
                else
                    errorbar(Dirs(idir), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                dir_resp(i,idir) = mean(resp(i,ind)-base(i,ind),2);
                hold on
            end
            if sum(h_dir(i,:),2)>0
                temp_resp = dir_resp(i,:);
                temp_resp(find(h_dir(i,:)==0)) = NaN;
                [max_val max_ind] = max(temp_resp,[],2);
                max_dir(i,:) = max_ind;
            else
                [max_val max_ind] = max(dir_resp(i,:),[],2);
                max_dir(i,:) = max_ind;
            end
            title([num2str(Dirs(max_dir(i,:))) ' deg'])
        end
    end
%         print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirTuning.pdf']),'-dpdf')

    
    nori = length(Dirs)/2;
    Ori = Dir;
    for iori = 1:nori
        ind = find(Dir == Dirs(iori+nori));
        Ori(ind) = Dirs(iori);
    end
    Oris = unique(Ori);
    h_ori = zeros(nCells, nori);
    p_ori = zeros(nCells, nori);
    ori_resp = zeros(nCells,nori);
    max_ori = zeros(nCells,1);
    figure;
    for i = 1:nCells
        if ~sum(no_match==i)
            subplot(n, n2, i)
            for iori = 1:nori
                if nOn>29
                    ind = find(Ori == Oris(iori));
                else
                    ind = find(Ori(1:ntrials-1) == Oris(iori));
                end
                [h_ori(i,iori), p_ori(i,iori)] = ttest(resp(i,ind), base(i,ind),'tail','right','alpha', 0.05/(nori-1));
                if h_ori(i,iori)
                    errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'or')
                else
                    errorbar(Oris(iori), mean(resp(i,ind)-base(i,ind),2),std(resp(i,ind)-base(i,ind),[],2)./sqrt(length(ind)),'ok')
                end
                ori_resp(i,iori) = mean(resp(i,ind)-base(i,ind),2);
                hold on
            end
            if sum(h_ori(i,:),2)>0
                temp_resp = ori_resp(i,:);
                temp_resp(find(h_ori(i,:)==0)) = NaN;
                [max_val max_ind] = max(temp_resp,[],2);
                max_ori(i,:) = max_ind;
            else
                [max_val, max_ind] = max(ori_resp(i,:),[],2);
                max_ori(i,:) = max_ind;
            end
            title([num2str(Oris(max_ori(i,:))) ' deg'])
        end
    end
    
    good_ind = unique([find(x)'; find(sum(h_dir,2)>0); find(sum(h_ori,2)>0)]);
%     print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuning.pdf']),'-dpdf')
    save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_trialData.mat']),'data_dfof','max_dir','h_dir', 'h_ori', 'max_ori','good_ind','dir_resp')
end

%% ori fitting
nOn = input.nScansOn;
nOff = input.nScansOff;
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
nTrials = length(dir_mat);
input.trialSinceReset = nTrials;

down = 10;
nframes = size(npSub_tc,1)./down;
nCells = size(npSub_tc,2);
data_tc_down = squeeze(mean(reshape(npSub_tc, [down,nframes,nCells]),1));
tuningDownSampFactor = down;

[avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC] = ...
    getOriTuningLG(data_tc_down,input,tuningDownSampFactor);
    vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:));

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningAndFits.mat']),...
            'avgResponseEaOri','semResponseEaOri','vonMisesFitAllCellsAllBoots','fitReliability','R_square', 'tuningTC')

%%
dir_mat = celleqel2mat_padded(input.tGratingDirectionDeg);
ori_mat = dir_mat;
ori_mat(find(dir_mat>=180)) = ori_mat(dir_mat>=180)-180;
oris = unique(ori_mat);
figure; 
if nCells<49
    [n n2] = subplotn(nCells);
else
    [n, n2] = subplotn(49);
end
start = 1;
x = 0;
for ic = 1:nCells
    if start > 49
        suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
        print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningFits_cells' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')
        start = 1;
        x = x+1;
        figure;
    end
    subplot(n,n2,ic-(x.*49))
    errorbar(oris,avgResponseEaOri(ic,:), semResponseEaOri(ic,:),'-o')
    hold on
    plot(0:180,vonMisesFitAllCellsAllBoots(:,1,ic));
    tit_str = num2str(chop(R_square(1,ic),2));
    if fitReliability(ic)<22.5
        tit_str = [tit_str '- R'];
    end
    title(tit_str)
    start = start+1;
end
suptitle([mouse ' ' date ' n = ' num2str(length(find(fitReliability<22.5))) '/' num2str(nCells) '- well-fit'])
print(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningFits' num2str(start-49) '-' num2str(start-1) '.pdf']),'-dpdf','-fillpage')

[max_resp prefOri] = max(vonMisesFitAllCellsAllBoots,[],1);
prefOri = squeeze(prefOri)-1;
prefOri_bootdiff = abs(prefOri(2:end,:)-prefOri(1,:));
prefOri_bootdiff(find(prefOri_bootdiff>90)) = 180-prefOri_bootdiff(find(prefOri_bootdiff>90));
ind_theta90 = find(prctile(prefOri_bootdiff,90,1)<22.5);
edges = [0 22.5:45:180]; 
[bin ind_bin] = histc(prefOri(1,:),edges);
ind_bin(find(ind_bin==5)) = 1;
bin(1) = bin(1)+bin(5);
bin(5) = [];

tunedCells = cell(1,length(bin));
for i = 1:length(bin)
    tunedCells{i} = intersect(find(ind_bin==i),ind_theta90);
end

save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_oriTuningInfo.mat']),...
    'prefOri', 'prefOri_bootdiff', 'ind_theta90', 'tunedCells');

%% identify cells that are red

% rename variables
fov_avg{1} = data_reg_avg;
dfmax{1} = data_dfof_max;
corrmap{1} = pix;
masks{1} = mask_cell;
corrmap_norm{1} = uint8((corrmap{1}./max(corrmap{1}(:))).*255);
brightnessScaleFactor = 0.5;
fov_norm{1} = uint8((fov_avg{1}./max(fov_avg{1}(:))).*255);
fov_norm{1}(fov_norm{1} > (brightnessScaleFactor*255)) = brightnessScaleFactor*255;

% process the red channel from a 1000 frame run
% selec the correct date = either "1" for D1 or "2" for D2/3
irun = 1;
WL = '920';
ImgFolder = strvcat('001');
run = catRunName(ImgFolder, nrun);
imgMatFile = [ImgFolder '_000_000.mat'];
CD = fullfile(gl_fn, [mouse '\' date '_' mouse '\' ImgFolder(irun,:)]);
cd(CD);
load(imgMatFile);
% fprintf(['Reading run ' num2str(irun) '- ' num2str(info.config.frames) ' frames \r\n'])
data_temp = sbxread(imgMatFile(1,1:11),0,info.config.frames);
mkdir(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run]));

if size(data_temp,1)>1
data_rg = squeeze(data_temp(1,:,:,:));
data_rr = squeeze(data_temp(2,:,:,:));
[out data_g_reg] = stackRegister(data_rg,mean(data_rg,3));
[out2 data_r_reg] = stackRegister_MA(data_rr,[],[],out);
red = mean(data_r_reg,3);
greenChImg = mean(data_g_reg,3);
clear data_temp 
fov_red{1} = uint8((red./max(red(:))).*255);
end

%% select red cells
% size of cell box
clear input
close all
w=30;
h=30;
buf = 3;

% get cell centroids
redGreenCells = struct;

cellPosition = regionprops(masks{1});
nc = length(cellPosition);

xCenter = cellfun(@(a) round(a(1)),{cellPosition.Centroid});
yCenter = cellfun(@(a) round(a(2)),{cellPosition.Centroid});

% index cells NOT too close to edge and NOT in black part of transformation
[ypix,xpix] = size(fov_avg{1});
goodCells = xCenter>(w/2) & xCenter<xpix-(w/2) & yCenter>(h/2) & yCenter<ypix-(h/2);

goodCells = goodCells & ...
    arrayfun(@(x) sum(sum(fov_norm{1}(masks{1}==x)))>0,1:nc);

    % fine register each cell    
mask_exp = zeros(size(fov_avg{1}));
mask_all = zeros(size(fov_avg{1}));

for icell = 1:nc
    if goodCells(icell)
        % find best shift
        day1_cell_avg = fov_norm{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        corr = corrmap_norm{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        day1_cell_max = dfmax{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        day1_red_avg = fov_red{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);
        
        mask = masks{1}(...
            yCenter(icell)-(h/2):yCenter(icell)+(h/2)-1,...
            xCenter(icell)-(w/2):xCenter(icell)+(w/2)-1);

            pass = true;
            figure;
            movegui('center')
            start = 1;
            subplot(1,4,start)
            imagesc(day1_cell_avg);axis image
            hold on
            bound = cell2mat(bwboundaries(mask(:,:,1)));
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);axis image
            title('avg')
            subplot(1,4,start+1)
            imagesc(corr);axis image
            hold on
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);axis image
            title('corr')
            subplot(1,4,start+2)
            imagesc(day1_cell_max);axis image
            hold on
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);
            title('dfof')
            subplot(1,4,start+3)
            imagesc(day1_red_avg);axis image
            hold on
            plot(bound(:,2),bound(:,1),'-','color','r','MarkerSize',2);
            title('Red')
            
            prompt = 'Choose one: 1- good, 2-okay, 3- none: ';
            x = input(prompt);
            switch x
                case 1
                    pass = true;
                    faint = false;
                case 2
                    pass = false;
                    faint = true;
                case 3
                    pass = false;
                    faint = false;
            end
    if pass
        redGreenCells(icell).pass = pass;
        redGreenCells(icell).faint = false;  
    elseif faint
        redGreenCells(icell).pass = false;
        redGreenCells(icell).faint = faint;  
    else 
        redGreenCells(icell).pass = false;
        redGreenCells(icell).faint = false;  
    end
    close all
    end  
end

goodCells = find([redGreenCells.pass]);
okayCells = find([redGreenCells.faint]);
redCells = sort([goodCells okayCells]);
redChImg = fov_red{1};
save(fullfile(fnout, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_multiday_alignment.mat']),'goodCells','okayCells','redCells','redChImg');
