% same thing as twoP_running_spk_graphs, !!! but using result of spike
% identification by doing deconvolution
% 2 photon: RUNNING ---- spikes
%plot average spike rate for all cells across trials, aligned with running onset and offset
%scatter plot running vs. still (each dot should be a cell)

%% assign document paths and experimental sessions
clear;
sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'}; 
days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
color_code = {'b','r','k','c'};

%% SECTION IV FR RUNNING VS. STAYTIONARY
for i = 4
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
    behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    threshold = -4;% the threshold you used in deconvolveCa function
    spk_deconv_output = load([image_analysis_dest sessions{i},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
    spk_inx_cl = spk_deconv_output.spk_inx_cl; % each cell is the indexes of CS for a PC
    FRstay_cell_cl = spk_deconv_output.FRstay_cell_cl;
    frm_run_cell = behav_output.frames_run_cell;
    frm_run = cell2mat(frm_run_cell);
    spkrun = zeros(1,size(spk_inx_cl,2));
    for n = 1: size(spk_inx_cl,2) %for each cell
        %count how many spikes are in the running periods: get the intersection
        %of indices of peaks above the threshold and indices of frames for running
        spkrun(n) = length(intersect(spk_inx_cl{n},frm_run));
    end
    t_run = length(frm_run)/30;
    FRrun_cells = spkrun/t_run;
    aveFR_run = mean(FRrun_cells);
    aveFR_stay = mean(FRstay_cell_cl);
    
    eg_stayvsRun = figure;
    scatter(FRstay_cell_cl,FRrun_cells,15,'filled','MarkerEdgeColor',[0.45 0.45 0.45],'MarkerFaceColor',[0.45 0.45 0.45]); hold on;
    scatter(aveFR_stay, aveFR_run,15,'filled','MarkerEdgeColor',[0.8431 0.0980 0.1098],'MarkerFaceColor',[0.8431 0.0980 0.1098]);hold on;
    xlabel('stationary'); ylabel('running');
    xlim([0,1.6]); ylim([0,1.6]);
    line = refline(1,0);
    line.Color = 'r';
    line.LineWidth = 1;
    axis square;
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',7);
    eg_stayvsRun.Units = 'centimeters';
    eg_stayvsRun.Position = [1 1 4 4];
    fig_name = 'eg_scatter_stayVsrun_190603_img1025';
    path = 'Z:\Analysis\figures\figure3_2Prun_spike\';
    orient(eg_stayvsRun,'landscape');
    print(eg_stayvsRun,[path,fig_name],'-r600','-depsc');
    
    
%     title([sessions{i} 'mean firing rate per second']);
%     savefig([image_analysis_dest sessions{i} '_meanFRstayVSrun_deconvolve' num2str(threshold)]);
%     save([image_analysis_dest sessions{i},'_spk_deconvolve_threshold' num2str(threshold) '.mat'],...
%         'FRrun_cells','aveFR_stay','aveFR_run','-append');
    
    %pause;
end

%finish
%% SECTION V : run trigger average and run offset average--- spikes rates
for i = 4
    image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
    behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
    frames_run_cell = behav_output.frames_run_cell;
    speed = double(behav_output.speed);
    frms_runTrig_mat = behav_output.frms_runTrig_mat;
    frms_runoff_mat = behav_output.frms_runoff_mat;
    avespd_runtrig_100ms = behav_output.avespd_runtrig_100ms;
    avespd_runoff_100ms = behav_output.avespd_runoff_100ms;
    stespd_runtrig_100ms = behav_output.stespd_runtrig_100ms;
    stespd_runoff_100ms = behav_output.stespd_runoff_100ms;
    
    threshold = -4;
    spk_deconv_output = load([image_analysis_dest sessions{i},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
    spk_bi_cellmat = spk_deconv_output.spk_logic_cl;
    spk_inx_cl = spk_deconv_output.spk_inx_cl;

    %run this script after dfOvF plot is done.
    
%     period = 9; % # of frames before and after running
%     befoRunStay = 30; %1s, # of frames that speed =0 before running
%     totalT = 60; %2s # of frmaes before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
%     aftRunOff = 30; %1s,# of frames that speed = 0 after running
%     [~,~,frms_runTrig_mat,frms_runoff_mat,~]= findFrames_runWindows_2P(speed,...
%         frames_run_cell,period,befoRunStay,totalT,aftRunOff);
    
    
    % aligned with running onset --------------------------------------------------------------------------------------------------
    %get a matrix of 0s and 1s, frames*windows*cells
    spk_runtrig_mat = zeros(size(frms_runTrig_mat,1),size(frms_runTrig_mat,2),size(spk_inx_cl,2));
    for n = 1: size(spk_inx_cl,2)                                                    %for each cell
        for j = 1:size(frms_runTrig_mat,2)                                    %for each running window
            spk_runtrig_mat (:,j,n) = spk_bi_cellmat(frms_runTrig_mat(:,j),n);%spk_runtrig_mat: frame*window*cell
        end
    end
    %calculate mean and ste
    spk_runtrig_cells = squeeze(mean(spk_runtrig_mat,2))*30;                    % average across windows, *30: average firing rate = average spike probablity * sampling rate
    spk_runtrig_grand = squeeze(mean(spk_runtrig_cells,2));                     % average across cells
    ste_runtrig = std(spk_runtrig_cells,0,2)/sqrt(size(spk_runtrig_cells,2));
%     % speed
%     speed_runtrig = speed(frms_runTrig_mat);
%     aveSpd_runtrig = mean(speed_runtrig,2);
%     steSpd_runtrig = std(speed_runtrig,0,2)/sqrt(size(speed_runtrig,2));
    %plot
    if size(frms_runTrig_mat,2)<=1
        disp('!!There is only one or no running trig window in fullfills the criteria, not plotting!!');
        disp(sessions{i});
    else
        x = (1:size(frms_runTrig_mat,1))/30-1;
        eg_spk_runtrig = figure; subplot(2,1,1);
        %shadedErrorBar(x,spk_runtrig_grand,ste_runtrig,{'color',[0.4660, 0.6740, 0.1880]}); % matlab default green
        shadedErrorBar(x,spk_runtrig_grand,ste_runtrig,{'color',[0.1922 0.6392 0.3294]});hold on;
        ylim([-0.1,2.6]); xlim([-1.1 1.6]);hold on;
        vline(0,'k'); ylabel('firing rate(Hz)'); 
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontSize',7);
        %title([sessions{i} ' aligned to running onset-deconvolution']);
        subplot(2,1,2);%plot(speed_runtrig);
        shadedErrorBar(x,avespd_runtrig_100ms*2*3.1415926*7.5/128, stespd_runtrig_100ms*2*3.1415926*7.5/128,{'color','k'});hold on;
        ylabel('speed(cm/s)'); xlabel('time from running onset(s)');
        ylim([0 15]);xlim([-1.1 1.6]);
        vline(0,'k'); 
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontSize',7);
        eg_spk_runtrig.Units = 'centimeters';
        eg_spk_runtrig.Position = [1 1 5.5 5];
        fig_name = 'eg_spk_runtrig_190603_img1025';
        path = 'Z:\Analysis\figures\figure3_2Prun_spike\';
        orient(eg_spk_runtrig,'landscape');
        print(eg_spk_runtrig,[path,fig_name],'-r600','-depsc');
        %title([sessions ' running triggered average speed']);
        %savefig([image_analysis_dest sessions{i} '_runTrigAve_FRwSpeed_deconvolve' num2str(threshold)]);
    end
    
    % aligned with running offset ------------------------------------------------------------------------------
    spk_runoff_mat = zeros(size(frms_runoff_mat,1),size(frms_runoff_mat,2),size(spk_inx_cl,2));
    for n = 1: size(spk_inx_cl,2)                                                    %for each cell
        for j = 1:size(frms_runoff_mat,2)                                    %for each running window
            spk_runoff_mat (:,j,n) = spk_bi_cellmat(frms_runoff_mat(:,j),n);%spk_runtrig_mat: frame*window*cell
        end
    end
    spk_runoff_cells = squeeze(mean(spk_runoff_mat,2))*30;                    % average across windows, *30: firing rate
    spk_runoff_grand = squeeze(mean(spk_runoff_cells,2));                     % average across cells
    ste_runoff = std(spk_runoff_cells,0,2)/sqrt(size(spk_runoff_cells,2));
    
%     % speed
%     speed_runoff = speed(frms_runoff_mat);
%     aveSpd_runoff = mean(speed_runoff,2);
%     steSpd_runoff = std(speed_runoff,0,2)/sqrt(size(speed_runoff,2));
    
    %plot
    if size(frms_runoff_mat,2)<=1
        disp('!!There is only one or no running off window in fullfills the criteria, not plotting!!');
        disp(sessions{i});
    else
        x = (1:size(frms_runoff_mat,1))/30-1;
        eg_spk_runoff = figure; subplot(2,1,1);
        shadedErrorBar(x,spk_runoff_grand,ste_runtrig,{'color',[0.1922 0.6392 0.3294]});hold on;
        ylim([-0.1,2.6]); xlim([-1.1 1.6]);hold on;
        vline(0,'k'); ylabel('firing rate(Hz)'); 
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontSize',7);
        %title([sessions{i} ' aligned to running offset-deconvolution']);
        subplot(2,1,2);
        shadedErrorBar(x,avespd_runoff_100ms*2*3.1415926*7.5/128, stespd_runoff_100ms*2*3.1415926*7.5/128,{'color','k'});hold on;
        ylabel('speed(cm/s)'); xlabel('time from running offset(s)');
        ylim([0 15]);xlim([-1.1 1.6]);
        vline(0,'k'); 
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'FontSize',7);
        eg_spk_runoff.Units = 'centimeters';
        eg_spk_runoff.Position = [1 1 5.5 5];
        fig_name = 'eg_spk_runoff_190603_img1025';
        path = 'Z:\Analysis\figures\figure3_2Prun_spike\';
        orient(eg_spk_runoff,'landscape');
        print(eg_spk_runoff,[path,fig_name],'-r600','-depsc');
        %savefig([image_analysis_dest sessions{i} '_runOffAve_FRwSpeed' num2str(threshold)]);
    end
    %pause;
%     save([image_analysis_dest sessions{i},'_spk_deconvolve_threshold' num2str(threshold) '.mat'],...
%         'spk_runtrig_mat','spk_runoff_mat','-append');
end


%%
% %% SECTION V : run trigger average and run offset average--- spikes rates
% % REGARDLESS OF whether the animal is stationary all the time before/after running
% for i = 1: size(sessions,2)
%     image_analysis_dest = [image_analysis_base, sessions{i}, '\'];
%     behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{i} '\'];
%     behav_output = load([behav_dest days{i} '_behavAnalysis.mat']);
%     frames_run_cell = behav_output.frames_run_cell;
%     speed = double(behav_output.speed);
%     threshold = -4;
%     spk_deconv_output = load([image_analysis_dest sessions{i},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
%     spk_bi_cellmat = spk_deconv_output.spk_logic_cl;
%     spk_inx_cl = spk_deconv_output.spk_inx_cl;
%     
%     period = 9; % # of frames before and after running
%     befoRunStay = 30; %1s, # of frames that speed =0 before running
%     totalT = 60; %2s # of frmaes before running onset + # of frames following running onset/ # of frames before running offset + # of frames after running offset
%     aftRunOff = 30; %1s,# of frames that speed = 0 after running
%     [~,~,frms_runTrig_NC_mat,frms_runoff_NC_mat,~]= findFrames_run_no_criteria_2P(speed,...
%         frames_run_cell,period,befoRunStay,totalT,aftRunOff);
%     % aligned with running onset ----------------------------------------------------------------------------
%     %get a matrix of 0s and 1s, frames*windows*cells
%     spk_runtrig_NC_mat = zeros(size(frms_runTrig_NC_mat,1),size(frms_runTrig_NC_mat,2),size(spk_inx_cl,2));
%     for n = 1: size(spk_inx_cl,2)                                                    %for each cell
%         for j = 1:size(frms_runTrig_NC_mat,2)                                    %for each running window
%             spk_runtrig_NC_mat (:,j,n) = spk_bi_cellmat(frms_runTrig_NC_mat(:,j),n);%spk_runtrig_mat: frame*window*cell
%         end
%     end
%     %calculate mean and ste
%     spk_runtrig_NC_cells = squeeze(mean(spk_runtrig_NC_mat,2))*30;                    % average across windows, *30: firing rate
%     spk_runtrig_NC_grand = squeeze(mean(spk_runtrig_NC_cells,2));                     % average across cells
%     ste_runtrig_NC = std(spk_runtrig_NC_cells,0,2)/sqrt(size(spk_runtrig_NC_cells,2));
%     % speed
%     speed_runtrig_NC = speed(frms_runTrig_NC_mat);
%     aveSpd_runtrig_NC = mean(speed_runtrig_NC,2);
%     steSpd_runtrig_NC = std(speed_runtrig_NC,0,2)/sqrt(size(speed_runtrig_NC,2));
%     %plot
%     x = (1:size(frms_runTrig_NC_mat,1));
%     figure; subplot(2,1,1);
%     shadedErrorBar(x,spk_runtrig_NC_grand,ste_runtrig_NC);hold on;
%     %ylim([-0.1,2]); hold on;
%     vline(31,'r'); ylabel('average firing rate'); xlabel('frame');
%     title([sessions{i} ' aligned to running onset INCLUSIVE-deconvolve']);
%     subplot(2,1,2);%plot(speed_runtrig);
%     shadedErrorBar(x,aveSpd_runtrig_NC, steSpd_runtrig_NC);hold on;
%     vline(31,'r'); ylabel('speed'); xlabel('frame');
%     savefig([image_analysis_dest sessions{i} '_runTrigAve_FRwSpeed_INCLUSIVE_deconvolve' num2str(threshold)]);
%     
%     
%     % aligned with running offset ------------------------------------------------------------------------------
%     spk_runoff_NC_mat = zeros(size(frms_runoff_NC_mat,1),size(frms_runoff_NC_mat,2),size(spk_inx_cl,2));
%     for n = 1: size(spk_inx_cl,2)                                                    %for each cell
%         for j = 1:size(frms_runoff_NC_mat,2)                                    %for each running window
%             spk_runoff_NC_mat (:,j,n) = spk_bi_cellmat(frms_runoff_NC_mat(:,j),n);%spk_runtrig_mat: frame*window*cell
%         end
%     end
%     spk_runoff_NC_cells = squeeze(mean(spk_runoff_NC_mat,2))*30;                    % average across windows, *30: firing rate
%     spk_runoff_NC_grand = squeeze(mean(spk_runoff_NC_cells,2));                     % average across cells
%     ste_runoff_NC = std(spk_runoff_NC_cells,0,2)/sqrt(size(spk_runoff_NC_cells,2));
%     
%     % speed
%     speed_runoff_NC = speed(frms_runoff_NC_mat);
%     aveSpd_runoff_NC = mean(speed_runoff_NC,2);
%     steSpd_runoff_NC = std(speed_runoff_NC,0,2)/sqrt(size(speed_runoff_NC,2));
%     
%     %plot
%     x = (1:size(frms_runoff_NC_mat,1));
%     figure; subplot(2,1,1);
%     shadedErrorBar(x,spk_runoff_NC_grand,ste_runoff_NC);hold on;
%     %ylim([-0.1,3]); hold on;
%     vline(31,'r'); ylabel('average firing rate');xlabel('frame');
%     title([sessions{i} ' aligned to running offset INCLUSIEVE-deconvolve']);
%     
%     subplot(2,1,2);%plot(speed_runoff);
%     shadedErrorBar(x,aveSpd_runoff_NC*2*3.1415926*7.5/128, steSpd_runoff_NC*2*3.1415926*7.5/128);hold on;
%     vline(31,'r'); ylabel('speed(cm/s)');xlabel('frame');%title([sessions ' running triggered average speed']);
%     
%     savefig([image_analysis_dest sessions{i} '_runOffAve_FRwSpeed_INCLUSIVE_deconvolve' num2str(threshold)]);
%     
%     save([image_analysis_dest sessions{i},'_spk_deconvolve_threshold' num2str(threshold) '.mat'],...
%         'spk_runtrig_NC_mat','speed_runtrig_NC','spk_runoff_NC_mat','speed_runoff_NC','-append');
% 
% 
% end

%% -------------------------------------------NOT USEFUL--------------plot all of the first derivatives: not very useful ------------------------------------------------------------
% %make a for loop to plot the shit so that the lines are not on top of each other
% colord=[         0         0    1.0000
%     0    0.4000         0
%     1.0000         0         0
%     0    0.7500    0.7500
%     0.7500         0    0.7500
%     0.8, 0.5, 0
%     0         0    0.5
%     0         0.85      0];
% 
% shift  = 0;
% fig = figure;
% for j = 1:size(deriv,2)
%     plot(deriv(:,j)+shift,'Color',colord(mod(j-1,size(colord,1))+1,:)); hold on
%     %     if ismember(i, [5,10,15,20,25,30,35,40,45,50])
%     %         hline(shift);
%     %     end
%     shift = shift+10000;  %10000 for 2P
% end
% set(gca,'YTick',10000:10000:shift);
% set(gca,'YTicklabel',(1:size(deriv,2)));
