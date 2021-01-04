%Analyze raw 2P data
%motion correction; PCA, ICA
%smooth and threshold mask, seperate overlaps and combine the one have a high correlation
%get TCs and plot masks

%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = '200320_img1063_tone'; 
%ID = '1016';
image_source_base  = 'Z:\Data\2photon\'; %location of permanently stored image files for retreiving meta data
%image_analysis_base    = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\'; 
image_analysis_base    = 'Z:\Analysis\motorizedWheel_Analysis\tone\imaging_analysis\'; 
%image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; 
%image_analysis_base    = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\'; 
%image_analysis_base    = 'Z:\Analysis\motorizedWheel_Analysis\airpuff\imaging_analysis\'; 

%image_source = [image_source_base, sessions,'\',ID,'\'];
image_source = [image_source_base, sessions,'\'];
image_analysis_dest = [image_analysis_base, sessions, '\' 'getTC\'];
%% 
cd(image_source);
order = '000';
%file = [ID '_000_' order];
file = [sessions '_000_' order];
%% Motion correct
%First, find a stack of frames that do not shift, this is the refrence for motion correct.
% read part of the sbx file, write into tiff, this allows you to look at
% the data in imageJ. Look at the tiff stack, find a stack of about 100
% frames that doesn't shift by eyeballing. 

% write tiff to get an idea about what the data looks like --------------------------------------------------------
%the very first frame in sbxfile should be 0, but the first frame of the
%new 2P is half black, so read it from the second frame
%!!!!!!!make sure you delete the speed of the first frame in behavior
frame_start = 1; 
nframes = 29999;%total frame number - frame_start
imgread = squeeze(sbxread(file,frame_start,nframes));
f1 = 2;
f2 = 1000;
img_tiff = imgread(:,:,f1:f2);
writetiff(img_tiff,[image_analysis_dest sessions '_'...
    order '_tiff_', num2str(f1), '_', num2str(f2)]);
%saveastiff(img_tiff,[image_analysis_dest sessions '_' order '_tiff_', num2str(f1) '_' num2str(f2)]);

%average every 500 frames and make figure, pick ref frame -----------------------------------------------------------
[npw,nph,nt] = size(imgread);
fprintf('Finding reference frame for motion correction...');
ncat = mod(-nt,500);
data = cat(3,imgread,nan(npw,nph,ncat));
img_refstacks = nanmean(permute(reshape(data,npw,nph,500,[]), [1 2 4 3]),4);
figure;
for i = 1:60
    subplot(8,8,i);imagesc(img_refstacks(:,:,i)); colormap gray;
    set(gca,'xticklabel',[],'yticklabel',[]);
    text(.8,.1,num2str(i),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end
savefig([image_analysis_dest sessions '_' order 'ref_frames']);
ref_frame = input('Input most stable region for reference: ');
img_ref = img_refstacks(:,:,ref_frame);
save([image_analysis_dest sessions '_' order '_img_rgs_ref' num2str(ref_frame) '.mat'], 'ref_frame','-v7.3');
%fprintf('done \n');

%choose stable frms by eyeball:
%stable_int = (180:330);img_ref = mean(imgread(:,:,stable_int),3);
%writetiff(img_ref,[image_analysis_base,sessions,'\' sessions '_'...
 %   order '_avetiffref_', num2str(stable_int(1)), '_', num2str(stable_int(end))]);
%motion register (use stack register, needs images that need to be registered 
% and reference of what these images need to be aligned to. Output: how
% much each frame moved and registed images ---------------------------------------------------
[rgs_out,img_rgs] = stackRegister(imgread, img_ref);
% look at the part of the registered image and save the new array ----------------------------------------------
gap = 50;
idx = (1:gap:nframes);
writetiff(img_rgs(:,:,idx),[image_analysis_dest sessions '_'...
    order '_rgstr_tiff_', num2str(frame_start), '_', num2str(nframes) '_',num2str(gap) '_ref' num2str(ref_frame)]);

save([image_analysis_dest sessions '_' order '_img_rgs_ref' num2str(ref_frame) '.mat'], 'img_rgs', 'rgs_out', 'img_ref','-append');
%another way to look at the registered movies: wirte an index(e.g.: all running onsets), and draw 100
%frames around each element, average each of those 100 frames and write a movie

%% PCA: reduce the dimensions that we don't need (the dimensions that do not
%contribute much to the varience of the data)
%first, down sample the raw data by averaging every 3 frames (reduces the
%temporal resolution but Jake said it's helpful for the PCA analysis

%img_ave = stackGroupProject(img_rgstr, 3); %average together every 3 frames
%imgsize = size(img_ave,3);
%img_rgs = img_rgs(:,:,101:end);
imgsize = size (img_rgs,3);
nPCA = 200; % the command window tells you "First #PCs contain #% of variance
            % if the nPCA is too good (contains all of the variances, ICA
            % will give you warinings (matrix is cloase to singular or
            % badly saled) and the # of convergence in rounds will be very
            % low like 3.
%run the PCA
[mixedsig_PCA, mixedfilters_PCA, CovEvals_PCA, ~, ~, ~] = CellsortPCA_2P_Jin(img_rgs,[1 imgsize], nPCA,[], []);
% if don't down sample, first 300PCs contain 34.3% of the variance. If down
% sample, first 300PCs contain 91.7% of the variance. Nathan said this number doesn't matter
% visualize 
figure; plot(CovEvals_PCA); % looks like can take 2-60th PCA or sth
savefig([image_analysis_dest sessions '_' order '_PCA', num2str(nPCA) 'CoEvals']);

figure;
for i = 1:196
subplot(14,14,i); imagesc(mixedfilters_PCA(:,:,i));
set(gca,'xticklabel',[],'yticklabel',[]);
    text(.8,.1,num2str(i),'fontsize',12,'color','w','fontweight','bold','unit','norm');
colormap gray
end
savefig([image_analysis_dest sessions '_' order '_PCA', num2str(nPCA)]);

save([image_analysis_dest sessions '_' order '_PCA_variables_', num2str(nPCA),'.mat'], 'mixedsig_PCA', 'mixedfilters_PCA', 'CovEvals_PCA', 'nPCA');


%% ICA: seperates independent spatial and temporal components
%PCuse =       1:125;
PCuse =       1:size(mixedfilters_PCA,3);
mu =          0.3; % weight of temporal info in spatio-teporal ICA
nIC =         200; % cannot be bigger than nPCA. If CoEvals doesn't change in later ICs, it will not converge!
ica_A_guess = []; %If this is empty than matlab will randomdize it and you can get different results, can see the random number generator in CellsortICA2P
termtol =      1e-6;
maxrounds =   3000;
npw =         264;
nph =         796;
[ica_sig, mixedfilters_ICA, ica_A, numiter] = CellsortICA_2P(mixedsig_PCA,...
                                                             mixedfilters_PCA,...
                                                             CovEvals_PCA, ...
                                                             PCuse,...
                                                             mu,...
                                                             nIC,...
                                                             [],...
                                                             termtol,...
                                                             maxrounds);
                                     
 icasig = permute(mixedfilters_ICA,[2,3,1]);% make it the same dimension as the raw data (pixel by pixel by time)
 ICA_variables.mu = mu; ICA_variables.nIC = nIC;  ICA_variables.termtol = termtol; ICA_variables.naxrounds = maxrounds;
 ICA_variables.npw = npw;  ICA_variables.nph = nph;
 save([image_analysis_dest sessions '_' order '_nIC' num2str(nIC) '_variables_for_nPCA_', ...
     num2str(nPCA), '.mat'], 'ica_sig', 'icasig', 'ICA_variables');
 writetiff(icasig, [image_analysis_dest sessions '_'...
    order '_icasig_for_nPCA', num2str(nPCA) '_nIC' num2str(nIC)]);
 figure; imagesc(sum(icasig,3));
 savefig([image_analysis_dest sessions '_' order '_icasig_sum_', num2str(nPCA) '_nIC' num2str(nIC)]);

 
%select which ICs to keep based on morphology
%disp(['Beginning IC selection for ', sessions, ' ',order, ' nPCA=', num2str(nPCA)])
%IC_use = IC_manual_check(icasig);
%save([out_dir, 'ICA_variables_for_nPCA_', num2str(nPCA), '.mat'], 'IC_use', '-append');

%% smooth and threshold the mask
%stack filter - acts as a low pass spatial filter, reduces noise. It will apply the Gaussian filter, and filter out low spatial frequency noise such as small blips and bloops in the ICs which do not belong to dendrites.
icasig_filt = stackFilter(icasig);

%set threshold a threshold for which pixels to include in a given dendrite's mask.
nIC = size(icasig_filt, 3);
cluster_threshold = 96.5; % this is using the top 3 percent of the fluorescence values, so brightest 3% is yes (1), and the rest is no (0)
%tried lower values for the threshold and turns out to have some wierd
%masks
mask_cell = zeros(size(icasig_filt));

for ic = 1:nIC
    %initialize sm_logical
    sm_logical = zeros(npw,nph);
    %convert to a binary mask (0 and 1)
    icasig_filt(:,:,ic) = imclearborder(icasig_filt(:,:,ic));
    sm_logical((icasig_filt(:,:,ic)> mean([max(prctile(icasig_filt(:,:,ic),cluster_threshold,1)) max(prctile(icasig_filt(:,:,ic),cluster_threshold,2))])))=1;
    sm_logical((icasig_filt(:,:,ic)<=mean([max(prctile(icasig_filt(:,:,ic),cluster_threshold,1)) max(prctile(icasig_filt(:,:,ic),cluster_threshold,2))])))=0;
    sm_logical = logical(sm_logical);
    %bwlabel identifies unconnected objects within a single mask. So 0=background 1=Object#1 2=Object#2
    if sum(sum(sm_logical,2),1) <51
        sm_logical = 0;
    end
    mask_cell(:,:,ic) = bwlabel(sm_logical); 
end
save([image_analysis_dest sessions '_' order, '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_mask.mat'], 'mask_cell');

%visualize the masks
figure('rend', 'painters', 'pos', [50 150 (796*1.5) (264*1.5)]); imagesc(sum(mask_cell,3));...
    title([sessions, order, '_', ' nPCA ', num2str(nPCA), ' mu ', num2str(mu), ' nIC ', num2str(nIC)]);
savefig([image_analysis_dest sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_cell_sum.fig']);

figure;
for i = 1:100
subplot(10,10,i); imagesc(mask_cell(:,:,i));
colormap gray
end
savefig([image_analysis_dest sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_cell_indi.fig']);


%% separate overlapping masks and combine the ones that have a high correlation
%split individual masks, remove small masks, deal with overlapping
mask_final = processMask(mask_cell);
mask_raw = reshape(mask_final, npw, nph);
figure; imagesc(mask_raw); truesize;
savefig([image_analysis_dest sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_process.fig']);

% combine highly correlated ICs to one
threshold = 0.8; %got the same thing when threshold = 0.8 and 0.9, tried different values, doesn't seem to change things
[ ~, mask3D, ~] = finalMask_Jin(img_rgs, mask_final, threshold);
%figure; imshow([image_analysis_dest 'AVG_' sessions '_' order '_rgstr_tiff_0_' num2str(nframes) '_50_jpeg.jpg']); hold on;
figure; imshow([image_analysis_dest 'AVG_' sessions '_' order '_rgstr_tiff_1_' num2str(nframes) '_50_ref' num2str(ref_frame) '_jpeg.jpg']); hold on;
for i  = 1:size(mask3D,3)
    bound = cell2mat(bwboundaries(mask3D(:,:,i)));
    randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
            num2str(i), 'color', 'y', 'FontSize', 8);
end
%figure; imagesc(sum(mask3D,3)); truesize; % got the same thing as the figure above
savefig([image_analysis_dest sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), ...
    '_coor' num2str(threshold),'_mask_wdendrites.fig']);
save([image_analysis_dest sessions '_' order, '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_coor' num2str(threshold) '_mask3D.mat'], 'mask3D');


%% look at the masks, mananully delete the ones that don't look like dendrites
mask3D(:,:,[]) = [];
%figure; imshow([image_analysis_dest 'AVG_' sessions '_' order '_rgstr_tiff_0_' num2str(nframes) '_50_jpeg.jpg']); hold on;
figure; imshow([image_analysis_dest 'AVG_' sessions '_' order '_rgstr_tiff_1_' num2str(nframes) '_50_ref' num2str(ref_frame) '_jpeg.jpg']); hold on;
for i  = 1:size(mask3D,3)
    bound = cell2mat(bwboundaries(mask3D(:,:,i)));
    randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), ...
            num2str(i), 'color', 'y', 'FontSize', 8);
end
hold on;
title([sessions ' masks_final']);
savefig([image_analysis_dest sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), ...
    '_coor' num2str(threshold),'_mask_wdendrites_final.fig']);
save([image_analysis_dest sessions '_' order, '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_coor' num2str(threshold) '_mask3D_final.mat'], 'mask3D');


%% get TCs
nmask = size(mask3D,3);
FrameRate = 30;
tc_avg = getTC(img_rgs, mask3D, nmask);

% delete the first frame !!!!!! make sure you also delete the first frame
% of speed. doing this because the first frame is half-black

% tc_avg = tc_avg((2:end),:);

rgstr_sum = sum(img_rgs,3);
plotTC_Jin(tc_avg, mask3D, rgstr_sum, 1:size(tc_avg,2), FrameRate);
% if plotting all the neurons make all lines look flat, plot some of them
% plotTC_Jin(tc_avg,~,~,1:40,~);
savefig([image_analysis_dest sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_TC.fig']);
mask_flat = plotMask_Jin(mask3D,1); %second input is either 0 or 1, if 1, makes a figure. 
savefig([image_analysis_dest sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_mask_plotMask.fig']);

data_corr = corrcoef(tc_avg);
figure; fig = imagesc(data_corr);
savefig([image_analysis_dest sessions '_' order, '_nPCA', num2str(nPCA),...
    '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold), '_corrcoef.fig']);

save([image_analysis_dest sessions '_' order, '_nPCA', ...
    num2str(nPCA), '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', ...
    num2str(cluster_threshold), '_TCave.mat'], 'tc_avg');

% %plot TC zoomed in if have too many neurons on the TC plot
% PC1 = 32;
% PC2 = 65;
% plotTC_Jin(tc_avg,mask3D, rgstr_sum,PC1:PC2,FrameRate); hold on;
% title([sessions ' TCave']);
% savefig([image_analysis_dest sessions '_' order, '_nPCA', num2str(nPCA),...
%     '_mu', num2str(mu), '_nIC', num2str(nIC), '_thresh', num2str(cluster_threshold),'_' num2str(PC1) '-' num2str(PC2) '_TC.fig']);
% 

%plotTC_Jin(tc_avg,1,1,32:65,1); this also works

