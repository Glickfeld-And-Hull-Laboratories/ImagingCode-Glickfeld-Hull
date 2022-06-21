clc; clear all; close all;
doRedChannel = 0;
ds = 'RFMapping_15Hz_ExptList_SG';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
%%
for iexp = 1
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).rfFolder;
time = expt(iexp).rfTime;
nrun = length(ImgFolder);
% run_str = catRunName(cell2mat(ImgFolder), nrun);
run_str = 'runs-002';

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';

fprintf(['2P imaging eye analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load data

load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))
load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))

%% eyetracking
calib = 1/26.6; %mm per pixel

% Load and combine eye tracking data
nrun = length(ImgFolder);
data = [];
for irun =  1:nrun
    CD = [SG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
    cd(CD);
    fn = [ImgFolder{irun} '_000_000_eye.mat'];

    data_temp = load(fn);          % should be a '*_eye.mat' file
    data_temp = squeeze(data_temp.data);

    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time{irun} '.mat'];
    load(fName);
    nFrames = input.counterValues{end}(end);
    data = cat(3, data, data_temp(:,:,1:nFrames));      % the raw images...
end

figure;
data_avg = mean(data,3);
imagesc(data_avg);
movegui('center')
ax = gca;
rect = getrect(ax);
datat = data_avg(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
figure;
imagesc(datat)
movegui('center')
while 1  % till broken out of

    % interactively get clicks
    [X Y selectionType] = getAPoint(gca);

    if isnan(X)
        key = lower(Y);
        switch key
          case char(13) % return
            break;  % out of loop, done
          case 'z' 
            imagesc(datat)
            rect = getrect(ax);
            datat = data(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3));
            imagesc(datat)
        end
        continue
    end
end
close all
data = data(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:);

%%
rad_range = [5 20];
warning off;
A = cell(size(data,3),1);
B = cell(size(data,3),1);
C = cell(size(data,3),1);
D = cell(size(data,3),1);
for n = 1:size(data,3)
    A{n} = [0,0];
    B{n} = [0];
    C{n} = [0];
    D{n} = [0];
end
eye = struct('Centroid',A,'Area',B,'Val',C,'SNR',D);
radii = [];
for n = 1:size(data,3)
    [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range,'Sensitivity',0.95);
              % pick the circle with best score
    if(isempty(center))
        eye(n).Centroid = [NaN NaN];    % could not find anything...
        eye(n).Area = NaN;
        eye(n).Val = NaN;
        eye(n).SNR = NaN;
    else
        snr = zeros(1,size(center,1));
        for idx = 1:size(center,1)
            t = double(data(:,:,n));
            vector_of_y_values = (1:size(data,1)) - center(idx,2);
            vector_of_x_values = (1:size(data,2)) - center(idx,1);
            [Yg, Xg] = ndgrid(vector_of_y_values, vector_of_x_values);
            idx1 = find(Xg.^2 + Yg.^2 < (radii(idx)/2).^2);
            idx2 = find(Xg.^2 + Yg.^2 < (radii(idx).*2.5).^2 & Xg.^2 + Yg.^2 > (radii(idx).*1.5).^2);
            snr(idx) = mean(t(idx1))./mean(t(idx2));
        end
        [v,idx] = max(snr);
        val = metric(idx);
        t = double(data(:,:,n));
        vector_of_y_values = (1:size(data,1)) - center(idx,2);
        vector_of_x_values = (1:size(data,2)) - center(idx,1);
        [Yg, Xg] = ndgrid(vector_of_y_values, vector_of_x_values);
        idx1 = find(Xg.^2 + Yg.^2 < (radii(idx)/2).^2);
        idx2 = find(Xg.^2 + Yg.^2 < (radii(idx).*2.5).^2 & Xg.^2 + Yg.^2 > (radii(idx).*1.5).^2);
        snr = mean(t(idx1))./mean(t(idx2));
        eye(n).SNR = snr;
        eye(n).Val = val;
        eye(n).Centroid = center(idx,:);
        eye(n).Area = pi*radii(idx)^2;
    end
    if mod(n,100)==0
        fprintf('Frame %d/%d\n',n,size(data,3));
    end
end
Centroid = cell2mat({eye.Centroid}');
Area = cell2mat({eye.Area}');
Val = double(cell2mat({eye.Val}'));
SNR = double(cell2mat({eye.SNR}'));
Eye_data = data;
nanframes(1,iexp) = length(find(isnan(Area)));

% no measurement frames
figure; 
subplot(2,2,1)
hist(sqrt(Area./pi));
xlabel('radius')
subplot(2,2,2)
hist(SNR);
xlabel('SNR')
subplot(2,2,3)
hist(Val);
xlabel('Metric')
movegui('center')

x1 = find(isnan(Area));
x2 = find(~isnan(Area));
x3 = unique([find(Val<0.1); find(Val<0.20 & SNR<1.7)]);

x = unique([x1; x3]);
if length(x)>25
    minx = 25;
else
    minx = length(x);
end

frames = sort(randsample(length(x),minx));
figure;
start = 1;
for i = 1:minx
    subplot(5,5,start);
    imagesq(data(:,:,x(frames(i)))); 
    hold on;
    scatter(Centroid(x(frames(i)),1), Centroid(x(frames(i)),2))
    title([num2str(chop(SNR(x(frames(i))),2)) ' ' num2str(chop(Val(x(frames(i))),2))])
    %title(num2str(x(frames(i))))
    start = start+1;
end
movegui('center')
sgtitle(['No pupil detected- ' num2str(length(x)) ' frames'])
print(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_noPupil2.pdf']),'-dpdf','-fillpage');

x = setdiff(x2,x3);
if length(x)>25
    minx = 25;
else
    minx = length(x);
end
frames = sort(randsample(length(x),minx));
figure;
start = 1;
for i = 1:minx
    subplot(5,5,start);
    imagesq(data(:,:,x(frames(i)))); 
    hold on;
    scatter(Centroid(x(frames(i)),1), Centroid(x(frames(i)),2))
    title([num2str(chop(SNR(x(frames(i))),2)) ' ' num2str(chop(Val(x(frames(i))),2))])
    %title(num2str(x(frames(i))))
    start = start+1;
end
movegui('center')
sgtitle('Pupil detected')
print(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Pupil.pdf']),'-dpdf','-fillpage');
        

%%
    
    %align eyetracking to 
     %reset frame counter    
     if exist('input.cStimOneOn', 'var')
        cStimOn = celleqel2mat_padded(input.cStimOneOn);
     else
         counter = celleqel2mat_padded(input.counter);
         cStimOn = counter + 1 - double(nOn+nOff); %magic number
     end 

     if ~exist('prewin_frames', 'var')   
        prewin_frames = nOff;
        postwin_frames = nOn;
     else
     end
        
    nanrun = ceil(500*(frame_rate/1000));
    Rad_temp = sqrt(Area./pi);
    Centroid_temp = Centroid;
    Rad_temp(unique([x1; x3]),:) =nan(length(unique([x1; x3])),1);
    Centroid_temp(unique([x1; x3]),:) = nan(length(unique([x1; x3])),2);
    sz = size(Eye_data);
    rad_mat_start = zeros(prewin_frames+postwin_frames, nTrials);
    centroid_mat_start = zeros(prewin_frames+postwin_frames,2, nTrials);
    eye_mat_start = zeros(sz(1), sz(2), prewin_frames+postwin_frames, nTrials);
   
    nframes = size(Rad_temp,1);
    
    for itrial = 1:nTrials
        if itrial == nTrials
            crange = [double(cStimOn(itrial))-prewin_frames:nframes];
        else
            crange = [double(cStimOn(itrial))-prewin_frames: double(cStimOn(itrial+1)-prewin_frames-1)];
        end
        if sum(isnan(Rad_temp(crange,1)),1)>0
            if sum(isnan(Rad_temp(crange,1)),1)./length(crange)> 0.25
                Rad_temp(crange,1) = NaN(length(crange),1);
                Centroid_temp(crange,:) = NaN(length(crange),2);
            else
                nanind = intersect(crange,find(isnan(Rad_temp)));
                dataind = intersect(crange,find(~isnan(Rad_temp)));
                for inan = 1:length(nanind)
                    gap = min(abs(nanind(inan)-dataind),[],1);
                    good_ind_stim = find(abs(nanind(inan)-dataind) == gap);
                    Rad_temp(nanind(inan),1) = mean(Rad_temp(dataind(good_ind_stim),1),1);
                    Centroid_temp(nanind(inan),:) = mean(Centroid(dataind(good_ind_stim),:),1);
                end
            end
        end
        if itrial < nTrials
            rad_mat_start(:,itrial) = Rad_temp(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
            centroid_mat_start(:,:,itrial) = Centroid_temp(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
            eye_mat_start(:,:,:,itrial) = Eye_data(:,:,cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1);
        else
            if (cStimOn(itrial)+postwin_frames)<nframes
                rad_mat_start(:,itrial) = Rad_temp(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
                centroid_mat_start(:,:,itrial) = Centroid_temp(cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1,:);
                eye_mat_start(:,:,:,itrial) = Eye_data(:,:,cStimOn(itrial)-prewin_frames:cStimOn(itrial)+postwin_frames-1);
            else
                rad_mat_start(:,itrial) = nan(prewin_frames+postwin_frames,1);
                centroid_mat_start(:,:,itrial) = nan(prewin_frames+postwin_frames,2,1);
                eye_mat_start(:,:,:,itrial) = nan(sz(1),sz(2),prewin_frames+postwin_frames,1);
            end
        end
            
    end
    rad_mat_calib = bsxfun(@times, rad_mat_start, calib);
    centroid_mat_calib = bsxfun(@times,centroid_mat_start,calib);
    t = mean(centroid_mat_calib(prewin_frames+1:end,:,:),1);
    rad_base = mean(rad_mat_calib(1:prewin_frames,:),1);
    rad_stim = mean(rad_mat_calib(prewin_frames+1:end,:),1);
    centroid_base = squeeze(mean(centroid_mat_calib(1:prewin_frames,:,:),1))./0.025;
    centroid_stim = squeeze(mean(centroid_mat_calib(prewin_frames+1:end,:,:),1))./0.025;

    figure; subplot(2,1,1)
    scatter(centroid_stim(1,:),centroid_stim(2,:), [], rad_stim); colorbar
    ind = find(~isnan(centroid_stim(1,:)));
    %centroid_med = geometric_median(centroid_stim(:,ind));
    centroid_med = findMaxNeighbors(centroid_stim(:,ind),2);
    hold on;
    scatter(centroid_med(1),centroid_med(2),'og')
    centroid_dist = sqrt((centroid_stim(1,:)-centroid_med(1)).^2 + (centroid_stim(2,:)-centroid_med(2)).^2);
    title('Color- radius')
    xlabel('x-pos')
    ylabel('y-pos')
    subplot(2,1,2)
    hist(centroid_dist,0:0.5:60)
    title([num2str(sum(centroid_dist<4)) ' trials w/in 4 deg'])
    sgtitle([num2str(sum(~isnan(centroid_dist))) '/' num2str(nTrials) ' measurable trials'])
    xlabel('Centroid distance from median')
    print(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupilPosDist.pdf']),'-dpdf','-fillpage');
    movegui('center')
        
    [n edges bin] = histcounts(centroid_dist,[0:2:30]);
    
    i = find(n);
    [n1 n2] = subplotn(length(i)); 
    figure;
    for ii = 1:length(i)
        subplot(n1,n2,ii)
        ind = find(bin== i(ii),1);
        if ii == 1
            ind_i = ind;
        end
        imagesc(mean(eye_mat_start(:,:,prewin_frames+1:end,ind),3))
        hold on
        plot(squeeze(nanmean(centroid_mat_start(prewin_frames+1:end,1,ind),1)), squeeze(nanmean(centroid_mat_start(prewin_frames+1:end,2,ind),1)),'or')
        plot(squeeze(nanmean(centroid_mat_start(prewin_frames+1:end,1,ind_i),1)), squeeze(nanmean(centroid_mat_start(prewin_frames+1:end,2,ind_i),1)),'ok')
        title([num2str(edges(ii)) '- ' num2str(length(find(bin== i(ii))))])
    end
    sgtitle('Example eye image by distance from median')
    print(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupilImgByDist.pdf']),'-dpdf','-fillpage');
    movegui('center')
   
    save(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']), 'rect', 'Area', 'Centroid', 'SNR', 'Val', 'frame_rate' , 'rad_mat_start','centroid_mat_start', 'cStimOn', 'rad_base','rad_stim','centroid_base', 'centroid_stim', 'centroid_dist', 'centroid_med');

end