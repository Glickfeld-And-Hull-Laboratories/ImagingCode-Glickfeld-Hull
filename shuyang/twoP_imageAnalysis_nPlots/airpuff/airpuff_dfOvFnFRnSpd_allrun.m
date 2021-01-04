%This script is to plot the PC activity and speed before and after airpuff onset. 
%In order to be included in frms_airstim_stay and frms_airstim_running,
%the bahevioral state can't change before and after the stimulus. 
%(i.e, if mice was running before airpuff, it has to be running as well after airpuff)
% this script will also create the variables that are useful:
% frms_airstim_run/frms_airstim_stay: trails*frames
% spk_airpuff_stay/run_mat: frames*trials*cells. 
% need a 3D matrix because need to keep the cell by cell variation. the average is the same when do average across trials first then across cells, or average cells first then across trials. However, the std is different. 
% FR_airstimAve_run: average firing rate across cells
% FR_airstim_AveSte_run: cell by cell variation.

%% Section I: set paths and create analysis folders for each session
%define the directory and files
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'191115_img1041'};
days = {'1041-191115_1'};
%sessionID = {'1023-190923_1','','','','','','','',''};
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
color_code = {'c','r','y','g'};

%%  sectionII: airpuff triggered response-spikes
for ii = 1: length(sessions)
    %% generate frame matrix based on the criteria, calculate average speed across trials
    behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    %behav_struct = load([behav_dest '\' days{ii} '_first18000frames_behavAnalysis.mat']);

    speed = behav_struct.speed; speed = double(speed);
    frames_run_cell = behav_struct.frames_run_cell;
    frames_stay_cell = behav_struct.frames_stay_cell;
    airpuffon = behav_struct.airpuffon1;
    %airpuffon = airpuffon(airpuffon<18000);
    
    period = 15; %#of frames before and after airpuff delivery
    
    % generate frame matrixes:how does neural activity change with airpuff delivery when the behvaior is the same?
    % frms_airstim_stay: 15 frames (0.5s)stationary before and after stim onset(the mice has to be stationary all the time!)
    % frms_airstim_allrun: !!9!! frames (0.43s)running before stim onset
    
    frms_airstim_stay = [];
    frms_airstim_allrun = [];
    for i = 1:length(airpuffon) 
        if airpuffon(i) - period < 0 || airpuffon(i)+period-1 > length(speed)
            continue
        elseif sum(speed((airpuffon(i)-period):(airpuffon(i)-1))==0) >= period && sum(speed(airpuffon(i):(airpuffon(i)+period-1))==0) >= period
            frms_airstim_stay = cat(1,frms_airstim_stay, airpuffon(i)- period : airpuffon(i)+period);% trail*frames
        elseif sum(speed((airpuffon(i)-period):(airpuffon(i)-1))~= 0) >= period-6 
            frms_airstim_allrun = cat(1,frms_airstim_allrun, airpuffon(i)- period : airpuffon(i)+period);% trial * frames
        end
    end
    %     % if the thershold for airpuff is shorter than 500ms
    %     for i = 2:length(airpuffon)-1 %ignoring the first and last stimulus because need a if loop for both(criteria same as below) and I have a lot of trials, can't fit them in the line below because than the index is out of dimension
    %         if (airpuffon(i+1)-airpuffon(i))< period ||(airpuffon(i)-airpuffon(i-1))<period %if airpuff is delivered more frequently than every 15 frames, get rid of these trials. These shouldn't happen because the threhold is always longer than period, but just in case if the threhold was set to a low value accidently
    %             continue
    %         elseif sum(speed((airpuffon(i)-period):(airpuffon(i)-1))==0) >= period && sum(speed(airpuffon(i):(airpuffon(i)+period-1))==0) >= period
    %             frms_airstim_stay = cat(1,frms_airstim_stay, airpuffon(i)- period : airpuffon(i)+period);% trail*frames
    %         elseif sum(speed((airpuffon(i)-period):(airpuffon(i)-1))~= 0) >= period-3 && sum(speed(airpuffon(i):(airpuffon(i)+period-1))~= 0) >= period-3
    %             frms_airstim_run = cat(1,frms_airstim_run, airpuffon(i)- period : airpuffon(i)+period);% trial * frames
%         end
%     end

  
 % generate matrix for speeds   
    spd_airstim_allrun = speed(frms_airstim_allrun);
    spd_airstim_stay = speed(frms_airstim_stay);
% averages and stes  
    spd_airstimAve_allrun = mean(spd_airstim_allrun);
    spd_airstimAveSte_allrun = std(spd_airstim_allrun)/sqrt(size(spd_airstim_allrun,1));
    
    spd_airstimAve_stay = mean(spd_airstim_stay);
    spd_airstimAveSte_stay = std(spd_airstim_stay)/sqrt(size(spd_airstim_stay,1));
    
%     save([behav_dest '\' days{ii} '_behavAnalysis.mat' ],...
%         'frms_airstim_allrun','frms_airstim_stay','spd_airstim_allrun',...
%         'spd_airstim_stay','spd_airstimAve_allrun','spd_airstimAve_stay',...
%         'spd_airstimAveSte_allrun','spd_airstimAveSte_stay','-append');
    
    %% spikes -- firing rate
    threshold = -4;
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    spk_deconv_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
    spk_logic = spk_deconv_output.spk_logic_cl;
    %spk_logic_ave = mean(spk_logic,2);
    %FR_ave = spk_logic_ave * 30; %firing rate = spike probablity*imaging rate
    
    % generate matrix for spikes
    spk_airpuff_allrun_mat = zeros(size(frms_airstim_allrun,1),size(frms_airstim_allrun,2),size(spk_logic,2)); %trial*frame*cell
    for i = 1: size(spk_logic,2)                                                    %for each cell
        for j = 1:size(frms_airstim_allrun,2)                                    %for each running trial
            spk_airpuff_allrun_mat (:,j,i) = spk_logic(frms_airstim_allrun(:,j),i);%spk_airpuff_run_mat: trial*frame*cell
        end
    end
    % average across trials:
    spk_airpuff_allrun_cells = squeeze(mean(spk_airpuff_allrun_mat)); %frame*cells
    %average across cells:
    FR_airstimAve_allrun = mean(spk_airpuff_allrun_cells,2)*30; %firing rate = spike probablity*imaging rate
    FR_airstimAveSte_allrun = std(spk_airpuff_allrun_cells,0,2)*30/sqrt(size(spk_airpuff_allrun_cells,2)); % cell by cell variation
    
    % generate matrix for spikes
    spk_airpuff_stay_mat = zeros(size(frms_airstim_stay,1),size(frms_airstim_stay,2),size(spk_logic,2));
    for i = 1: size(spk_logic,2)                                                    %for each cell
        for j = 1:size(frms_airstim_stay,2)                                    %for each running trial
            spk_airpuff_stay_mat (:,j,i) = spk_logic(frms_airstim_stay(:,j),i);%spk_airpuff_run_mat: trial*frame*cell
        end
    end
    % average across trials:
    spk_airpuff_stay_cells = squeeze(mean(spk_airpuff_stay_mat));
    %average across cells:
    FR_airstimAve_stay = mean(spk_airpuff_stay_cells,2)*30; %firing rate = spike probablity*imaging rate
    FR_airstimAveSte_stay = std(spk_airpuff_stay_cells,0,2)*30/sqrt(size(spk_airpuff_stay_cells,2)); % cell by cell variation
    
%     save([image_analysis_dest '\' sessions{ii}, '_spk_deconvolve_threshold' num2str(threshold) '.mat' ],...
%         'spk_airpuff_allrun_mat','spk_airpuff_stay_mat','FR_airstimAve_allrun',...
%         'FR_airstimAveSte_allrun','FR_airstimAve_stay','FR_airstimAveSte_stay',...
%         'spk_airpuff_allrun_cells','spk_airpuff_stay_cells','-append');
 %{   
% plot only spike rate and speed
    % figure
    fig = figure;
    % if there's only one trial for running or stationary, plot that single
    % trial
    x = (1: length(spd_airstimAve_stay))/30;
    subplot(2,2,1);
    shadedErrorBar(x,FR_airstimAve_run,FR_airstimAveSte_run);
    title('airpuff stim running');  ylabel('firing rate');vline(period/30,'r');
    xlim([0,1]);ylim([0,4]);
    subplot(2,2,3);
    shadedErrorBar(x,spd_airstimAve_run*2*3.1415926*7.5/128,spd_airstimAveSte_run*2*3.1415926*7.5/128);
    xlabel('time(s)'); ylabel('speed(cm/s)'); vline(period/30,'r');
    xlim([0 1]);%ylim([0 20]);
    text(0.1,min(spd_airstimAve_run*2*3.1415926*7.5/128),['n = ' num2str(size(frms_airstim_run,1))]);
    %xlim([0,31]);
    
    subplot(2,2,2);
    shadedErrorBar(x,FR_airstimAve_stay,FR_airstimAveSte_stay);
    title('airpuff stim stay'); ylabel('firing rate');vline(period/30,'r');
    xlim([0,1]);ylim([0,4]);
    subplot(2,2,4);
    shadedErrorBar(x,spd_airstimAve_stay*2*3.1415926*7.5/128,spd_airstimAveSte_stay*2*3.1415926*7.5/128);
    xlabel('time(s)'); ylabel('speed(cm/s)'); vline(period/30,'r');
    xlim([0 1]);%ylim([-1 1])
    text(0.1,0.05,['n = ' num2str(size(frms_airstim_stay,1))]);
    %xlim([0,31]);
    
    supertitle(sessions{ii});
    
    saveas(fig, [image_analysis_dest '\' sessions{ii} '_airpuffTrigAve_FRnSpd']);    
  %}  
end



%%  sectionIII: airpuff triggered response-df/F
%frm_airstim_run is from the second section of this script, so have to run
%the section section first.
for ii = 1: length(sessions)
   % load data
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\']; 
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvF_btm = dfOvF_strct.dfOvF_btm_cl;
    %dfOvF_stay = dfOvF_strct.dfOvF_stay;
    behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    speed = behav_struct.speed; speed = double(speed);
    frms_airstim_stay = behav_struct.frms_airstim_stay;
    frms_airstim_allrun = behav_struct.frms_airstim_allrun;
    spd_airstimAve_allrun = behav_struct.spd_airstimAve_allrun;
    spd_airstimAveSte_allrun = behav_struct.spd_airstimAveSte_allrun;
    spd_airstimAve_stay = behav_struct.spd_airstimAve_stay;
    spd_airstimAveSte_stay = behav_struct.spd_airstimAveSte_stay;
 
 % use bottom 10% of fluorescence:
 % ------------------------------------------------------------------------------------------------------------------------------------------------------
    % generate matrix for dfOvF_btm during running
    dfOvFbtm_airpuff_allrun_mat = zeros(size(frms_airstim_allrun,1),size(frms_airstim_allrun,2),size(dfOvF_btm,2)); % trial*frame*cell
    for i = 1: size(dfOvF_btm,2)                                                    %for each cell
        for j = 1:size(frms_airstim_allrun,1)                                    %for each running trial
            dfOvFbtm_airpuff_allrun_mat(j,:,i) = dfOvF_btm(frms_airstim_allrun(j,:),i);%spk_airpuff_run_mat: trial*frame*cell
        end
    end
    % average across trials:
    dfOvFbtm_airpuff_allrun_cells = squeeze(mean(dfOvFbtm_airpuff_allrun_mat)); %frame*cells
    %average across cells:
    dfOvFbtm_airstimAve_allrun = mean(dfOvFbtm_airpuff_allrun_cells,2); %firing rate = spike probablity*imaging rate
    dfOvFbtm_airstimAveSte_allrun = std(dfOvFbtm_airpuff_allrun_cells,0,2)/sqrt(size(dfOvFbtm_airpuff_allrun_cells,2)); % cell by cell variation
    
     % generate matrix for dfOvF_btm during stationary
    dfOvFbtm_airpuff_stay_mat = zeros(size(frms_airstim_stay,1),size(frms_airstim_stay,2),size(dfOvF_btm,2));
    for i = 1: size(dfOvF_btm,2)                                                    %for each cell
        for j = 1:size(frms_airstim_stay,2)                                    %for each running trial
            dfOvFbtm_airpuff_stay_mat(:,j,i) = dfOvF_btm(frms_airstim_stay(:,j),i);%spk_airpuff_run_mat: trial*frame*cell
        end
    end
    % average across trials:
    dfOvFbtm_airpuff_stay_cells = squeeze(mean(dfOvFbtm_airpuff_stay_mat)); %frame*cells
    %average across cells:
    dfOvFbtm_airstimAve_stay = mean(dfOvFbtm_airpuff_stay_cells,2); 
    dfOvFbtm_airstimAveSte_stay = std(dfOvFbtm_airpuff_stay_cells,0,2)/sqrt(size(dfOvFbtm_airpuff_stay_cells,2)); % cell by cell variation
    
% ------------------------------------------------------------------------------------------------------------------------------------------------------  
  
% % use fluorescence during stationary without airpuff
%     % generate matrix for dfOvF_stay during running
%     dfOvFstay_airpuff_allrun_mat = zeros(size(frms_airstim_allrun,1),size(frms_airstim_allrun,2),size(dfOvF_stay,2));
%     for i = 1: size(dfOvF_stay,2)                                                    %for each cell
%         for j = 1:size(frms_airstim_allrun,2)                                           %for each running trial
%             dfOvFstay_airpuff_allrun_mat(:,j,i) = dfOvF_stay(frms_airstim_allrun(:,j),i);  %spk_airpuff_run_mat: trial*frame*cell
%         end
%     end
%     % average across trials:
%     dfOvFstay_airpuff_allrun_cells = squeeze(mean(dfOvFstay_airpuff_allrun_mat)); %frame*cells
%     %average across cells:
%     dfOvFstay_airstimAve_allrun = mean(dfOvFstay_airpuff_allrun_cells,2)*30; %firing rate = spike probablity*imaging rate
%     dfOvFstay_airstimAveSte_allrun = std(dfOvFstay_airpuff_allrun_cells,0,2)*30/sqrt(size(dfOvFstay_airpuff_allrun_cells,2)); % cell by cell variation
%     
%      % generate matrix for dfOvF_stay during stationary
%     dfOvFstay_airpuff_stay_mat = zeros(size(frms_airstim_stay,1),size(frms_airstim_stay,2),size(dfOvF_stay,2));
%     for i = 1: size(dfOvF_stay,2)                                                    %for each cell
%         for j = 1:size(frms_airstim_stay,2)                                    %for each running trial
%             dfOvFstay_airpuff_stay_mat(:,j,i) = dfOvF_stay(frms_airstim_stay(:,j),i);%spk_airpuff_run_mat: trial*frame*cell
%         end
%     end
%     % average across trials:
%     dfOvFstay_airpuff_stay_cells = squeeze(mean(dfOvFstay_airpuff_stay_mat)); %frame*cells
%     %average across cells:
%     dfOvFstay_airstimAve_stay = mean(dfOvFstay_airpuff_stay_cells,2)*30; %firing rate = spike probablity*imaging rate
%     dfOvFstay_airstimAveSte_stay = std(dfOvFstay_airpuff_stay_cells,0,2)*30/sqrt(size(dfOvFstay_airpuff_stay_cells,2)); % cell by cell variation
%     
%     save([image_analysis_dest '\' sessions{ii} '_dfOvF.mat' ],...
%         'dfOvFbtm_airpuff_allrun_mat','dfOvFbtm_airpuff_stay_mat',...
%         'dfOvFbtm_airpuff_allrun_cells','dfOvFbtm_airstimAve_allrun','dfOvFbtm_airstimAveSte_allrun',...
%         'dfOvFbtm_airpuff_stay_cells','dfOvFbtm_airstimAve_stay','dfOvFbtm_airstimAveSte_stay',...
%         '-append');
    % ------------------------------------------------------------------------------------------------------------------------------------------------------
   %{
    %plot
    period = 15; %#of frames before and after airpuff delivery
    fig = figure;
    x = (1:length(spd_airstimAve_stay))/30;
    subplot(3,2,1);
    shadedErrorBar(x,dfOvFbtm_airstimAve_run,dfOvFbtm_airstimAveSte_run);
    title('airpuff, running, df/F bottom 10%'); ylabel('df/f'); vline(period/30,'r');
    xlim([0,1]);
    subplot(3,2,3);
    shadedErrorBar(x,dfOvFstay_airstimAve_run,dfOvFstay_airstimAveSte_run);
    title('df/f F during stationary'); ylabel('df/f'); vline(period/30,'r');
    xlim([0,1]);
    subplot(3,2,5);
    shadedErrorBar(x,spd_airstimAve_run*2*3.1415926*7.5/128,spd_airstimAveSte_run*2*3.1415926*7.5/128);
    xlabel('time(s)'); ylabel('speed(cm/s)'); vline(period/30,'r');
    xlim([0 1]);
    text(0.1,min(spd_airstimAve_run*2*3.1415926*7.5/128),['n = ' num2str(size(frms_airstim_run,1))]);   
    
    subplot(3,2,2);
    shadedErrorBar(x,dfOvFbtm_airstimAve_stay,dfOvFbtm_airstimAveSte_stay);
    title('airpuff, stationary, df/F bottom 10%'); ylabel('df/f'); vline(period/30,'r');
    xlim([0,1]);
    subplot(3,2,4);
    shadedErrorBar(x,dfOvFstay_airstimAve_stay,dfOvFstay_airstimAveSte_stay);
    title('df/f F during stationary');ylabel('df/f'); vline(period/30,'r');
    xlim([0,1]);
    subplot(3,2,6);
    shadedErrorBar(x,spd_airstimAve_stay*2*3.1415926*7.5/128,spd_airstimAveSte_stay*2*3.1415926*7.5/128);
    xlabel('time(s)'); ylabel('speed(cm/s)'); vline(period/30,'r');
    xlim([0 1]);
    text(0.1,min(spd_airstimAve_stay*2*3.1415926*7.5/128),['n = ' num2str(size(frms_airstim_stay,1))]);
    
    supertitle(sessions{ii});
    saveas(fig, [image_analysis_dest '\' sessions{ii} '_airpuffTrigAve_dfOvFnSpd']);   
    %}
end

    %% plot
    % figure
    airpff_eg_fig = figure;
    % if there's only one trial for running or stationary, plot that single
    % trial
    
    x = (1: length(spd_airstimAve_stay))/30 - 0.5-(1/30);
    subplot(3,2,1);
    shadedErrorBar(x,dfOvFbtm_airstimAve_allrun,dfOvFbtm_airstimAveSte_allrun,{'color',[0.8431 0.0980 0.1098]});
    %title('airpuff stim running'); 
    ylabel('df/f');
    xlim([-0.5,0.5]);ylim([0,0.6]);
    vline((period+1)/30 - 0.5-(1/30),'k');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',8);
    subplot(3,2,3);
    shadedErrorBar(x,FR_airstimAve_allrun,FR_airstimAveSte_allrun,{'color',[0.1922 0.6392 0.3294]});
    ylabel('firing rate(Hz)');
    xlim([-0.5,0.5]); ylim([0 6]);
    vline((period+1)/30 - 0.5-(1/30),'k');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',8);
    subplot(3,2,5);
    shadedErrorBar(x,spd_airstimAve_allrun*2*3.1415926*7.5/128,spd_airstimAveSte_allrun*2*3.1415926*7.5/128,{'color','k'});
    xlabel('time from airpuff(s)'); ylabel('speed(cm/s)');
    %text(0.1,min(spd_airstimAve_allrun+10),['n = ' num2str(size(frms_airstim_allrun,1))]);
    xlim([-0.5,0.5]);ylim([-0.5 20]);
    vline((period+1)/30 - 0.5-(1/30),'k');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',8);
     
    subplot(3,2,2);
    shadedErrorBar(x,dfOvFbtm_airstimAve_stay,dfOvFbtm_airstimAveSte_stay,{'color',[0.8431 0.0980 0.1098]});
    %title('airpuff stim stay'); 
    %ylabel('df/f');
    xlim([-0.5,0.5]);ylim([0,0.6]);
    vline((period+1)/30 - 0.5-(1/30),'k');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',8);
    subplot(3,2,4);
    shadedErrorBar(x,FR_airstimAve_stay,FR_airstimAveSte_stay,{'color',[0.1922 0.6392 0.3294]});
    %ylabel('firing rate');
    xlim([-0.5,0.5]); ylim([0 6]);
    vline((period+1)/30 - 0.5-(1/30),'k');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',8);
    subplot(3,2,6);
    shadedErrorBar(x,spd_airstimAve_stay*2*3.1415926*7.5/128,spd_airstimAveSte_stay*2*3.1415926*7.5/128,{'color','k'});
     %ylabel('speed(cm/s)');
    %text(0.1,5,['n = ' num2str(size(frms_airstim_stay,1))]);
    xlim([-0.5,0.5]);ylim([-0.5 20]);
    vline((period+1)/30 - 0.5-(1/30),'k');
    xlabel('time from airpuff(s)');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',8);
    
    airpff_eg_fig.Units = 'centimeters';
    airpff_eg_fig.Position = [3 3 11 10];
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',8);
    fig_name = 'eg_191115_img1041_airpuff';
    path = 'Z:\Analysis\figures\figure4_airpuff_selfpace\';
    print(airpff_eg_fig,[path,fig_name],'-r600','-depsc');
    
    %suptitle(sessions{ii});
    
    %saveas(airpff_eg_fig, [image_analysis_dest '\' sessions{ii} '_airpuffTrigAve_allrun']);
    
