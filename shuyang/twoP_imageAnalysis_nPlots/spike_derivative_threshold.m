% spike identification using the standard derivations of the first
% derivative of raw traces. 
% function twoP_best_spkthreshold tells you what's the optimal std(i.e., 1 std ot 2 or 3), 
% then section II draws the GUI showing raw traces with identified events
% and first derivatives.
%% assign document paths and experimental sessions
clear;
% sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
% days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};
% image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; 

sessions = '190603_img1025'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; 
image_analysis_dest = [image_analysis_base, sessions,'\'];
image_analysis_dest_deriv = [image_analysis_base, sessions, '\derivative\'];
if ~exist(image_analysis_dest_deriv)
    mkdir(image_analysis_dest_deriv);
end

% behavior analysis results 
days = '1025-190603_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% SECTION I  SPIKE RATES
% get derivatives from raw fluorescence
TCave = load([image_analysis_dest sessions '_deconvolution_thresh-4_TCave_cl']); % use only good cells delected by deconvolution
TCave = TCave.TCave_cl;
deriv = diff(TCave,1,1); % derivatives shoud have positive and negative values
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_stay_cell = behav_output.frames_stay_cell;
frm_stay = cell2mat(frm_stay_cell);
%  decide the threshold for spikes,  get average number of spikes during stationary (bestFR)
[aveFRsWstds, best_thres,bestFR,std_best,spk_inx,FRstay_cells,...
    spk_bi_cellmat] = twoP_best_spkthreshold (deriv, frm_stay,TCave);
%savefig([image_analysis_dest_deriv sessions '_hist_stayFR_wdiffthresh']);
% 
% % plot spike logic
% figure;
% imagesc(spk_bi_cellmat);
% xlabel('Cell #');
% ylabel('Frame #');
% title (sessions);
% savefig([image_analysis_dest sessions '_scatter_spikeLogic']);
% 
% % plot spike rates for each cell during stationary
% figure; 
% scatter(1:length(FRstay_cells), FRstay_cells, 'bo');
% ylabel('Spike Rate (Hz)');
% xlabel('Cell #');
% title (sessions);
% savefig([image_analysis_dest sessions '_scatter_spikeRate']);

save([image_analysis_dest_deriv sessions,'_spikes.mat'],'aveFRsWstds','best_thres'...
    ,'bestFR','std_best','spk_inx','spk_bi_cellmat','FRstay_cells');

%% SECTION II plots to go back and look at the raw data

% % plot first derivatives with horizontal std lines --------------------------------------------------
% 
% listfrm = [1:300; 3301:3600;6801:7100;12051:12350;19001:19300;28501:28800; 16881:17180;  ...
%   21661:21960; 26001:26300 ];
% listfrm = listfrm';
% [fig_deriv] = GUI_w_hline(deriv, listfrm,std_deriv,std2,std2_5,std3);
% savefig([image_analysis_dest sessions '_GUIderiv_stdlines']);
% 
% % after function best_spkthreshold_2P, go back to raw fluorescence(TCave) and see if the time points of spikes make sense
% % plot red dots over rawTC when deriv > std_bestv for all cells and several frames---------------------------------------------------
% listfrm = [1:300; 3301:3600;6801:7100;12051:12350;19001:19300;28501:28800;   16881:17180;  ...
%   21661:21960; 26001:26300 ];
% listfrm = listfrm';
% [fig] = GUI_rawTrace_wSpkevents(TCave,listfrm,spk_inx,sessions);  
% savefig([image_analysis_dest sessions '_GUI_TCave_w_spkEvent_stdbest']); 
% 
% % make GUI with suplots: raw fluorescence with red circles and first derivatives
% % with std lines
% 
% [fig] = GUI_rawTrace_nderivThresh(TCave,deriv,spk_bi_cellmat,sessions);
% savefig([image_analysis_dest sessions '_GUI_TCave_w_deriv']); 
% 
% % Plot raw fluorescence for the running triggered average windows -----
% [fig_runtrig] = GUI_w_vline(TCave,frms_runTrig_mat,spk_inx); %!!!!! you will get an error if you have less than 9 running windows
% savefig([image_analysis_dest sessions '_GUIraw_fluo_RUNtrig']);
