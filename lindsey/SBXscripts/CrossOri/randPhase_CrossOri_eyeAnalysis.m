clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);
%%
for iexp = 9
mouse = expt(iexp).mouse;
date = expt(iexp).date;
area = expt(iexp).img_loc{1};
ImgFolder = expt(iexp).coFolder;
time = expt(iexp).coTime;
nrun = length(ImgFolder);
run_str = catRunName(cell2mat(ImgFolder), nrun);

LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

fprintf(['2P imaging eye analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
for irun=1:nrun
    fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
end

%% load data

load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))

%% eyetracking
calib = 1/26.6; %mm per pixel

% Load and combine eye tracking data
irun =  1;
CD = [LG_base '\Data\2P_images\' mouse '\' date '\' ImgFolder{irun}];
cd(CD);
fn = [ImgFolder{irun} '_000_000_eye.mat'];
data = load(fn);          % should be a '*_eye.mat' file

data = squeeze(data.data);      % the raw images...

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

rad_range = [3 19];
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
    [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range,'Sensitivity',0.9);
    [val,idx] = max(metric);          % pick the circle with best score
    if(isempty(center))
        eye(n).Centroid = [NaN NaN];    % could not find anything...
        eye(n).Area = NaN;
        eye(n).Val = NaN;
        eye(n).SNR = NaN;
    else
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

x1 = find(isnan(Area));
x2 = find(~isnan(Area));
x3 = find(Val<0.26 & SNR<1.9);

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
suptitle('No pupil detected')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_noPupil2.pdf']),'-dpdf','-fillpage');

x = setdiff(x2,x3);
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
suptitle('Pupil detected')
print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Pupil.pdf']),'-dpdf','-fillpage');
        

%%
    
    %align eyetracking to 
     %reset frame counter    
    cStimOn = celleqel2mat_padded(input.cStimOneOn);
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
    suptitle([num2str(sum(~isnan(centroid_dist))) '/' num2str(nTrials) ' measurable trials'])
    xlabel('Centroid distance from median')
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupilPosDist.pdf']),'-dpdf','-fillpage');

        
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
        title(num2str(edges(ii)))
    end
    suptitle('Example eye image by distance from median')
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupilImgByDist.pdf']),'-dpdf','-fillpage');

    
    centroid_dist_sf = cell(1,nSF);
    centroid_med_sf = cell(1,nSF);
    if nSF>1
        for isf = 1:nSF
            ind_plaid = intersect(find(maskCon_all == maskCons(2)), find(stimCon_all == stimCons(2)));
            ind_sf = intersect(ind_plaid,intersect(find(SF_all == SFs(isf)),find(~isnan(centroid_stim(1,:)))));
            centroid_med_sf{isf} = findMaxNeighbors(centroid_stim(:,ind),2);
            centroid_dist_sf{isf} = sqrt((centroid_stim(1,:)-centroid_med(1)).^2 + (centroid_stim(2,:)-centroid_med(2)).^2);
        end
    end    
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']), 'rect', 'Area', 'Centroid', 'SNR', 'Val', 'frame_rate' , 'rad_mat_start','centroid_mat_start', 'cStimOn', 'rad_base','rad_stim','centroid_base', 'centroid_stim', 'centroid_dist', 'centroid_med', 'centroid_dist_sf', 'centroid_med_sf' );
    close all
end
%     %% eye plots
% trial_n = zeros(nMaskCon,nStimCon,nMaskPhas);
% trialInd = cell(nMaskCon,nStimCon,nMaskPhas);
% for im = 1:nMaskCon
%     ind_mask = find(maskCon_all == maskCons(im));
%     for it = 1:nStimCon
%         ind_stim = find(stimCon_all == stimCons(it));
%         ind_sm = intersect(ind_mask,ind_stim);
%         if it>1 & im>1
%             for ip = 1:nMaskPhas
%                 ind_phase = find(maskPhas_all == maskPhas(ip));
%                 ind = intersect(ind_phase,ind_sm);
%                 trialInd{im,it,ip} = ind;
%                 trial_n(im,it,ip) = length(ind);
%             end
%         else
%             trialInd{im,it,1} = ind_sm;
%             trial_n(im,it,1) = length(ind);
%         end
%     end
% end
% 
% figure;
% start = 1;
% n = 1;
% for iC = 1:length(resp_ind)
%     iCell = resp_ind(iC);
%     if start > 25
%         suptitle({['All responsive cells- Test = ' num2str(stimCons(2)) ' Mask = ' num2str(maskCons(3))], ['Trials: ' num2str(length(find(~isnan(centroid_dist)))) ' measured; ' num2str(length(find(centroid_dist<5))) ' <5 deg; ' num2str(length(find(centroid_dist<2))) ' <2 deg']})
%         print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respByPhase_eyeTrack_allResp_' num2str(n) '.pdf']),'-dpdf','-fillpage');
%         start = 1;
%         figure;
%         n = n+1;
%     end
%     subplot(5,5,start)
%     test_avg = mean(resp_cell{1,2,1}(iCell,:),2);
%     mask_avg = mean(resp_cell{end,1,1}(iCell,:),2);
%     test_sem = mean(resp_cell{1,2,1}(iCell,:),2)./sqrt(size(resp_cell{1,end,1}(iCell,:),2));
%     mask_sem = mean(resp_cell{end,1,1}(iCell,:),2)./sqrt(size(resp_cell{end,1,1}(iCell,:),2));
%     shadedErrorBar(1:360,repmat(test_avg,[1 360]),repmat(test_sem,[1 360]));
%     hold on
%     shadedErrorBar(1:360,repmat(mask_avg,[1 360]),repmat(mask_sem,[1 360]));
%     resp_all = [];
%     stim_all = [];
%     for ip = 1:nMaskPhas
%         [memb ind] = ismember(trialInd{end,2,ip},find(~isnan(centroid_dist)));
%         resp_avg(1,ip) = mean(resp_cell{end,2,ip}(iCell,find(ind)),2);
%         resp_sem(1,ip) = std(resp_cell{end,2,ip}(iCell,find(ind)),[],2)./sqrt(size(resp_cell{end,2,ip}(iCell,find(ind)),2));
%         resp_all = [resp_all resp_cell{end,2,ip}(iCell,find(ind))];
%         stim_all = [stim_all ip.*ones(size(resp_cell{end,2,ip}(iCell,find(ind))))];
%     end
%     errorbar(maskPhas,resp_avg,resp_sem)
%     p1 = anova1(resp_all, stim_all,'off');
%     resp_all = [];
%     stim_all = [];
%     for ip = 1:nMaskPhas
%         [memb ind] = ismember(trialInd{end,2,ip},find(centroid_dist<5));
%         resp_avg(1,ip) = mean(resp_cell{end,2,ip}(iCell,find(ind)),2);
%         resp_sem(1,ip) = std(resp_cell{end,2,ip}(iCell,find(ind)),[],2)./sqrt(size(resp_cell{end,2,ip}(iCell,find(ind)),2));
%         resp_all = [resp_all resp_cell{end,2,ip}(iCell,find(ind))];
%         stim_all = [stim_all ip.*ones(size(resp_cell{end,2,ip}(iCell,find(ind))))];
%     end
%     errorbar(maskPhas,resp_avg,resp_sem)
%     p2 = anova1(resp_all, stim_all,'off');
%     resp_all = [];
%     stim_all = [];
%     for ip = 1:nMaskPhas
%         [memb ind] = ismember(trialInd{end,2,ip},find(centroid_dist<2));
%         resp_avg(1,ip) = mean(resp_cell{end,2,ip}(iCell,find(ind)),2);
%         resp_sem(1,ip) = std(resp_cell{end,2,ip}(iCell,find(ind)),[],2)./sqrt(size(resp_cell{end,2,ip}(iCell,find(ind)),2));
%         resp_all = [resp_all resp_cell{end,2,ip}(iCell,find(ind))];
%         stim_all = [stim_all ip.*ones(size(resp_cell{end,2,ip}(iCell,find(ind))))];
%     end
%     errorbar(maskPhas,resp_avg,resp_sem)
%     p3 = anova1(resp_all, stim_all,'off');
%     title([num2str(chop(p3,2))])
%     ylabel('dF/F')
%     xlabel('Phase (deg)')
%     start = start+1;
% end
% suptitle({['All responsive cells- Test = ' num2str(stimCons(2)) ' Mask = ' num2str(maskCons(3))], ['Trials: ' num2str(length(find(~isnan(centroid_dist)))) ' measured; ' num2str(length(find(centroid_dist<5))) ' <5 deg; ' num2str(length(find(centroid_dist<2))) ' <2 deg']})
% print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respByPhase_eyeTrack_allResp_' num2str(n) '.pdf']),'-dpdf','-fillpage');
% close all
% end