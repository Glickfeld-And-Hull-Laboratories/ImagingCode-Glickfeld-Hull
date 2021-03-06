%apply stackregister and PCA/ICA across sessions
file_info_CRP;
usFacs = 100;
currCD = cd;
frames_session_1 = num_good_frames(1);
frames_per_preview = 500;
%out_dir_base = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\combined_masks\';
out_dir_base = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\FoV_matching_outputs_thresholding\';

img_comb = [];
iii = 0;
for sub = [animal_to_be_analyzed]%size(mouseID,2)    %img93 day1 and post learning 500ms [22:25]   %img92 day1 and post learning [9,10]    %img91 day1/post 500ms [5:6]    %img90 day1/post 500ms [2:3]   %img89 day1/post 500ms [23,27]   %img94 day1/post500ms [31:32]
    iii = iii+1;
    for rID = 1:3;
        %session_runID = runID{rID};
        data_dir = fullfile('Z:\Data\2P_imaging',[dates{sub} '_' mouseID{sub}], mouseID{sub});
        config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
        img_fn   = dir(fullfile(data_dir,['*' runID{rID} '.sbx']));
        if ~isempty(img_fn) & exist([data_dir, '\', img_fn.name], 'file');
            load([data_dir, '\', config_fn.name]);
            nframes = info.config.frames;
            fprintf('loading sbx image file: %s %s\n', dates{sub}, img_fn.name);
            cd([data_dir, '\']);
            if num_good_frames(iii)>info.config.frames
                num_good_frames(iii)=info.config.frames;
                if iii==1
                    frames_session_1 = info.config.frames;
                end
            end
            %img = squeeze(sbxread(img_fn.name(1:end-4),0,num_good_frames(iii)));   %,num_good_frames(iii)
            [info, img, t] = load_sbx_movie(img_fn.name(1:end-4), data_dir, num_good_frames(iii));
            break
        else
            continue
        end
    end
    img_comb = cat(3, img_comb, img);
end
out_dir = [out_dir_base, mouseID{sub}, '\'];
if ~exist(out_dir);
    mkdir(out_dir);
end
clear img;
cd(currCD);

%% stackRegister
% find most stable reference image among the 100 randomly selected images
% from the first session. Then motion register both sessions to that image.
ref30 = img_comb(:,:,randi([1,frames_session_1],1,30)); %randomly selects frames to sample from session 1
[npw, nph, nt] = size(img_comb);
 
sf= round(linspace(1, frames_session_1, 100));
samp100 = img_comb(:,:,sf);
dshift = [];

for r = 1:size(ref30,3)
    [reg_out,aa] = stackRegister(samp100, ref30(:,:,r));
    dshift = [dshift;mean(((reg_out(:,3).^2)+(reg_out(:,4).^2)))];
end

min_f = find(dshift == min(dshift), 1, 'first');
img_ref = ref30(:,:,min_f);
[reg_out, img_reg] = stackRegister(img_comb, img_ref);
save([out_dir, 'Reg_out.mat'], 'reg_out','img_ref');

%write a preview tiff to compare motion reg, write max proj tiff as well
writetiff(img_comb(:,:,[1:500, frames_session_1+1:frames_session_1+501]), [out_dir, 'day1_postlearning_CoReg_500_per_session']); %--------------------------------------------------------
img_comb_max_proj = cat(3, max(img_comb(:,:,[1:500]),[],3), max(img_comb(:,:,[(frames_session_1+1):(frames_session_1+501)]),[],3));
writetiff(img_comb_max_proj, [out_dir, 'day1_postlearning_CoReg_Max_projects_for_each']);
clear img_comb

%% PCA/ICA

%compute principal components ---------------------------------
nPCA = 500; %100 for old datasets, 500 for newer
img_pca = img_reg(:,:,1:2:(frames_session_1/2)); % downsample in time by 2 or 5
nf = size(img_pca,3);
img_fn = img_fn.name;
[mixedsig, mixedfilters, CovEvals, ~, ~, ...
    ~] = CellsortPCA_2P(img_pca,[1 nf], nPCA,[], out_dir, img_fn, []);    %only use the first session to generate the mask
%             [mixedsig, mixedfilters, CovEvals, ~, ~, ~] = CellsortPCA2(img_reg,[1 nFrames], nPCA, [], out_dir, img_fn, []);

%compute independent components ------------------------------
PCuse = 1:nPCA;
mu = 0.98; % spatial temporal mix
if strcmp(mouseID{sub}, 'img90');
    nIC=400;
elseif strcmp(mouseID{sub}, 'img91');
    nIC=100; 
else
    nIC = 300;  %400- img90 100- img91
end
termtol = 0.00001;
maxrounds = 1000;
%ICA function
[ica_sig, mixedfilters, ica_A, numiter] = CellsortICA_2P(mixedsig, ...
    mixedfilters, CovEvals, PCuse, mu, nIC, [], termtol, maxrounds);

%save ICA outputs from session 1
icasig = permute(mixedfilters,[2,3,1]);
save([out_dir, 'ICA.mat'], 'ica_sig', 'icasig');

%% Create mask and process ICA outputs
nIC = size(icasig,3);
icasig = stackFilter(icasig);
mask_cell = zeros(size(icasig));
sm_logical = zeros(npw,nph);
if strcmp(mouseID{sub}, 'img93');
    cluster_threshold = 96;
else
    cluster_threshold = 95;
end

%thresholding to get each IC
for ic = 1:nIC
    icasig(:,:,ic) = imclearborder(icasig(:,:,ic));% remove ICs touching borders
    sm_logical((icasig(:,:,ic)>mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=1;
    sm_logical((icasig(:,:,ic)<=mean([max(prctile(icasig(:,:,ic),cluster_threshold,1)) max(prctile(icasig(:,:,ic),cluster_threshold,2))])))=0;
    sm_logical = logical(sm_logical);
    mask_cell(:,:,ic) = bwlabel(sm_logical);
end

% separate multiple ICs and remove small ones and deal with overlapping 
mask_final = processMask_CRP(mask_cell);
mask_raw = reshape(mask_final, npw, nph);
figure; imagesc(mask_raw); truesize
threshold = 0.8; % correlation threshold

%combine highly correlated ICs
[ ~, mask3D, ~] = finalMask(img_reg(:,:,1:10:end), mask_final, threshold, out_dir);
nmask = size(mask3D,3);

%% plot TCs
FrameRate = 30;
tc_avg = getTC(img_reg, mask3D, nmask);

saveData = 1;
reg_max = max(img_reg(:,:,[1:100]),[],3);
plotTC_CRP(tc_avg, mask3D, reg_max, 1:size(tc_avg,2), FrameRate, out_dir, saveData);
mask_flat = plotMask(mask3D, saveData, out_dir);

data_corr = corrcoef(tc_avg);
figure; fig = imagesc(data_corr);
%saveas(fig, [out_dir, 'data_corr.fig']);
%print([out_dir, 'data_corr.eps'],'-depsc')

mask_final = processMask(mask3D);
sz = [npw, nph];
save([out_dir, 'Results.mat'], 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr', 'sz');
%save([out_dir, 'ROI_TCs.mat'],'tc_avg', 'mask_final', 'sz');
save([out_dir, 'PCA_ICA_settings'], 'threshold', 'cluster_threshold', 'mu', 'nIC', 'termtol', 'maxrounds');
%close all

