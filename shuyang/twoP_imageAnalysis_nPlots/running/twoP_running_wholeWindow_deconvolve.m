%% 2 photon: RUNNING - 1s before running, whole running window, 1s after running (df/f and spikes)
% !! using result of spike identification by deconvolution.
% plots for df/f and firing rate with averages of each time points
% plots for df/f and firing data, topologize all trials.

%% assign document paths and experimental sessions
clear;
sessions = '190603_img1025'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1025-190603_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};

%% SECTION I df/f 
dfOvF = load([image_analysis_dest sessions '_dfOvF.mat']);
avedfOvF = dfOvF.avedfOvF;
stedfOvF = dfOvF.stedfOvF;
behav_output = load([behav_dest days '_behavAnalysis.mat']);
speed = double(behav_output.speed);
frames_run_cell = behav_output.frames_run_cell;
% get the frames of 1s before and 1s after all running trials
time = 30;%1s
frm_bfrun = zeros(size(frames_run_cell,2),time);
frm_aftrun = zeros(size(frames_run_cell,2),time);
for t = 1: size(frames_run_cell,2) %for each trial
    frm_bfrun(t,:)= (frames_run_cell{t}(1)-time:frames_run_cell{t}(1)-1);
    frm_aftrun(t,:) = (frames_run_cell{t}(end)+1:frames_run_cell{t}(end)+time);
end
% write frames_run_cell into a matrix, can't use findFrames_run_2P function
% because the matrix from that function only includes running trials that are longer than 1s
run_lengths = cell2mat(cellfun(@size,frames_run_cell, 'UniformOutput',0));
run_longest = max(run_lengths);
frms_run_alltials_mat = nan(size(frames_run_cell,2), run_longest);
for t = 1: size(frms_run_alltials_mat,1)
    temp = frames_run_cell{t};
    frms_run_alltials_mat(t,1:size(temp,2)) = temp;
end
%------------------------------------------------------------------------------------------------------------------------------------------
% put df/f into the matrix above and calculate mean and ste across trials.
% get aves and stes of different time windows and put them together
% can't put frames before run, run, and after run together first because then NaNs in the running window 
% will be at the end of each row, matlab doesn't automatically do integers, NaN, integers
dfOvF_runmat = nan(size(frms_run_alltials_mat,1),size(frms_run_alltials_mat,2));
for t = 1:size(frms_run_alltials_mat,1)
    temp = frms_run_alltials_mat(t,~isnan(frms_run_alltials_mat(t,:)));
    dfOvF_runmat(t,1:size(temp,2)) = avedfOvF(temp); 
end
mean_dfOvF_run = mean(dfOvF_runmat,'omitnan');
ste_dfOvF_run = std(dfOvF_runmat,0,1,'omitnan')/sqrt(size(dfOvF_runmat,1));

% calculate mean df/f before and after running
dfOvF_bfRun = avedfOvF(frm_bfrun);
mean_dfOvF_bfRun = mean(dfOvF_bfRun,'omitnan');
ste_dfOvF_bfRun = std(dfOvF_bfRun,0,1,'omitnan')/sqrt(size(dfOvF_bfRun,1));

dfOvF_aftRun = avedfOvF(frm_aftrun);
mean_dfOvF_aftRun = mean(dfOvF_aftRun,'omitnan');
ste_dfOvF_aftRun = std(dfOvF_aftRun,0,1,'omitnan')/sqrt(size(dfOvF_aftRun,1));

mean_dfOvF_wholeWindow = [mean_dfOvF_bfRun,mean_dfOvF_run,mean_dfOvF_aftRun];
ste_dfOvF_wholeWindow = [ste_dfOvF_bfRun,ste_dfOvF_run,ste_dfOvF_aftRun];

%----------------------------------------------------------------------------------------------------------------------------
%put speed into the running matrix
speed_runmat = nan(size(frms_run_alltials_mat,1),size(frms_run_alltials_mat,2));
for t = 1:size(frms_run_alltials_mat,1)
    temp = frms_run_alltials_mat(t,~isnan(frms_run_alltials_mat(t,:)));
    speed_runmat(t,1:size(temp,2)) = speed(temp); 
end
mean_speed_run = mean(speed_runmat,'omitnan');
ste_speed_run = std(speed_runmat,0,1,'omitnan')/sqrt(size(speed_runmat,1));

% calculate mean df/f before and after running
speed_bfRun = speed(frm_bfrun);
mean_speed_bfRun = mean(speed_bfRun,'omitnan');
ste_speed_bfRun = std(speed_bfRun,0,1,'omitnan')/sqrt(size(speed_bfRun,1));

speed_aftRun = speed(frm_aftrun);
mean_speed_aftRun = mean(speed_aftRun,'omitnan');
ste_speed_aftRun = std(speed_aftRun,0,1,'omitnan')/sqrt(size(speed_aftRun,1));

mean_speed_wholeWindow = [mean_speed_bfRun,mean_speed_run,mean_speed_aftRun];
ste_speed_wholeWindow = [ste_speed_bfRun,ste_speed_run,ste_speed_aftRun];

%------------------------------------------------------------------------------------------------------------------------------------
% plot
x = (1:length(mean_speed_wholeWindow))/30;
run_wholetime_dfOvF = figure;
subplot(2,1,1);
shadedErrorBar(x,mean_dfOvF_wholeWindow,ste_dfOvF_wholeWindow);
ylabel ('df/f');vline(time/30,'r'); 
vline((length(mean_dfOvF_wholeWindow)-time)/30,'r');
%xlim([0,length(mean_dfOvF_wholeWindow)+5]);
%xlim([7.33,10.33]);
%ylim([-0.3 0.1]);
%ylim([min(mean_dfOvF_wholeWindow)-0.05, max(mean_dfOvF_wholeWindow)+0.05]);
title(['df/f before, during, and after running' sessions]);

subplot(2,1,2);
shadedErrorBar(x,mean_speed_wholeWindow*2*3.1415926*7.5/128,ste_speed_wholeWindow*2*3.1415926*7.5/128);
ylabel ('speed(cm/s)');vline(time/30,'r'); 
vline((length(mean_dfOvF_wholeWindow)-time)/30,'r');
xlabel('time(s)');
%xlim([0,length(mean_speed_wholeWindow)+5]);
%xlim([7.33,10.33]); ylim([-50,50]);
%ylim([min(mean_speed_wholeWindow)-5,max(mean_speed_wholeWindow)+5]);
saveas(run_wholetime_dfOvF, [image_analysis_dest sessions, '_dfOvF_run_allWindows']);

save([image_analysis_dest sessions '_dfOvF.mat'],'mean_dfOvF_wholeWindow','ste_dfOvF_wholeWindow','mean_speed_wholeWindow','ste_speed_wholeWindow','-append');


%% firng rate 
%CAUTION: right now this session has to be run after the above session,
%because the code for generating the matrix of the frames and calculating mean speed are not re-written here.
threshold = -4;
spk_deconv_output = load([image_analysis_dest sessions,'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
spk_bi_cellmat = spk_deconv_output.spk_logic;
ave_spk_prob = mean(spk_bi_cellmat,2);
aveFR = ave_spk_prob*30; % average firing rate  = average spike probablity * sampling rate
aveFR = aveFR';

%---------------------------------------------------------------------------------------------------------------------------
%put aveFR into the running matrix
aveFR_runmat = nan(size(frms_run_alltials_mat,1),size(frms_run_alltials_mat,2));
for t = 1:size(frms_run_alltials_mat,1)
    temp = frms_run_alltials_mat(t,~isnan(frms_run_alltials_mat(t,:)));
    aveFR_runmat(t,1:size(temp,2)) = aveFR(temp); 
end
mean_aveFR_run = mean(aveFR_runmat,'omitnan');
ste_aveFR_run = std(aveFR_runmat,0,1,'omitnan')/sqrt(size(aveFR_runmat,1));

% calculate mean aveFR before and after running
aveFR_bfRun = aveFR(frm_bfrun);
mean_aveFR_bfRun = mean(aveFR_bfRun,'omitnan');
ste_aveFR_bfRun = std(aveFR_bfRun,0,1,'omitnan')/sqrt(size(aveFR_bfRun,1));

aveFR_aftRun = aveFR(frm_aftrun);
mean_aveFR_aftRun = mean(aveFR_aftRun,'omitnan');
ste_aveFR_aftRun = std(aveFR_aftRun,0,1,'omitnan')/sqrt(size(aveFR_aftRun,1));

mean_aveFR_wholeWindow = [mean_aveFR_bfRun,mean_aveFR_run,mean_aveFR_aftRun];
ste_aveFR_wholeWindow = [ste_aveFR_bfRun,ste_aveFR_run,ste_aveFR_aftRun];

%------------------------------------------------------------------------------------------------------------------------------------
% plot
x = (1: length(mean_aveFR_wholeWindow))/30;
run_wholetime_FR = figure;
subplot(2,1,1);
shadedErrorBar(x,mean_aveFR_wholeWindow,ste_aveFR_wholeWindow);
ylabel ('firing rate');vline(time/30,'r'); vline((length(mean_aveFR_wholeWindow)-time)/30,'r');
%xlim([0,length(mean_aveFR_wholeWindow)+5]);
%ylim([min(mean_aveFR_wholeWindow)-0.3,max(mean_aveFR_wholeWindow)+0.3]);
title(['firing rate before, during, and after running deconvolve' sessions]);

subplot(2,1,2);
shadedErrorBar(x,mean_speed_wholeWindow*2*3.1415926*7.5/128,ste_speed_wholeWindow*2*3.1415926*7.5/128);
ylabel ('speed(cm/s)');vline(time/30,'r'); vline((length(mean_aveFR_wholeWindow)-time)/30,'r');
xlabel('time(s)');
%xlim([0,length(mean_speed_wholeWindow)+5]);
%ylim([min(mean_speed_wholeWindow)-5,max(mean_speed_wholeWindow)+5]);
saveas(run_wholetime_FR, [image_analysis_dest sessions, '_FR_run_allWindows_deconvolve' num2str(threshold)]);

save([image_analysis_dest sessions '_spk_deconvolve_threshold' num2str(threshold) '.mat'],'mean_aveFR_wholeWindow','ste_aveFR_wholeWindow','mean_speed_wholeWindow','ste_speed_wholeWindow','-append');


%% topology - df/f 
% make df/f and speed into cell
dfOvF_run_cell = cell(size(frames_run_cell,1),size(frames_run_cell,2));
for c = 1: size(frames_run_cell,2)
    dfOvF_run_cell{c} = (avedfOvF(frames_run_cell{c}))';
end

speed_run_cell = cell(size(frames_run_cell,1),size(frames_run_cell,2));
for c = 1: size(frames_run_cell,2)
    speed_run_cell{c} = (speed(frames_run_cell{c}))';
end

% topology
resolution_ratio = 2000;   %divide the total process into how many parts
dfOvF_topo_mat=zeros(length(dfOvF_run_cell),resolution_ratio);  
for i=1:1:length(dfOvF_run_cell)
    dfOvF_temp = dfOvF_run_cell{i};
    mapping=linspace(1,length(dfOvF_temp)+1,resolution_ratio);
    for j=1:1:length(dfOvF_temp)
        dfOvF_topo_mat(i,(mapping>=j)) = dfOvF_temp(j);
    end 
end 
dfOvF_topo = mean(dfOvF_topo_mat);
ste_dfOvF_topo = std(dfOvF_topo_mat,0,1)/sqrt(size(dfOvF_topo_mat,1));


speed_topo_mat=zeros(length(speed_run_cell),resolution_ratio);  
for i=1:1:length(speed_run_cell)
    speed_temp=speed_run_cell{i};
    mapping=linspace(1,length(speed_temp)+1,resolution_ratio);
    for j=1:1:length(speed_temp)
        speed_topo_mat(i,(mapping>=j))=speed_temp(j);
    end 
end 
speed_topo = mean(speed_topo_mat);
ste_speed_topo = std(speed_topo_mat,0,1)/sqrt(size(speed_topo_mat,1));

dfOvF_topo_wholeWindow = [mean_dfOvF_bfRun, dfOvF_topo, mean_dfOvF_aftRun];
dfOvF_ste_topo_wholeWindow = [ste_dfOvF_bfRun,ste_dfOvF_topo, ste_dfOvF_aftRun];

speed_topo_wholeWindow = [mean_speed_bfRun, speed_topo, mean_speed_aftRun];
speed_ste_topo_wholeWindow = [ste_speed_bfRun, ste_speed_topo, ste_speed_aftRun];

% plot
x = [(1:1:time),linspace(1+time,time+length(mean_speed_run),resolution_ratio),(length(mean_speed_wholeWindow)-29:1:length(mean_speed_wholeWindow))];
x = x/30;
run_wholetime_dfOvF_topo = figure;
subplot(2,1,1);
shadedErrorBar(x,dfOvF_topo_wholeWindow,dfOvF_ste_topo_wholeWindow);
ylabel ('df/f');vline(time/30,'r'); vline((length(mean_speed_wholeWindow)-time)/30,'r');
%xlim([0,length(mean_speed_wholeWindow)+5]);
xlim([0 10.33]); 
ylim([-0.2 0.1]);
%ylim([min(dfOvF_topo_wholeWindow)-0.05,max(dfOvF_topo_wholeWindow)+0.05]);
title(['df/f topologized before, during, and after running' sessions]);

subplot(2,1,2);
shadedErrorBar(x,speed_topo_wholeWindow*2*3.1415926*7.5/128,speed_ste_topo_wholeWindow*2*3.1415926*7.5/128);
ylabel ('speed(cm/s)');vline(time/30,'r'); vline((length(mean_speed_wholeWindow)-time)/30,'r');
xlabel('time(s)');
%xlim([0,length(mean_speed_wholeWindow)+5]);
xlim([0 10.33]); 
ylim([-5 20]);
%ylim([min(speed_topo_wholeWindow)-5,max(speed_topo_wholeWindow)+5]);
saveas(run_wholetime_dfOvF_topo, [image_analysis_dest sessions, '_dfOvF_run_allWindows_topo']);

save([image_analysis_dest sessions '_dfOvF.mat'],'dfOvF_topo_wholeWindow','dfOvF_ste_topo_wholeWindow','speed_topo_wholeWindow','speed_ste_topo_wholeWindow','-append');


%% topology - firing rate
% make df/f and speed into cell
aveFR_run_cell = cell(size(frames_run_cell,1),size(frames_run_cell,2));
for c = 1: size(frames_run_cell,2)
    aveFR_run_cell{c} = (aveFR(frames_run_cell{c}))';
end

% topology
resolution_ratio = 2000;   %divide the total process into how many parts
aveFR_topo_mat=zeros(length(aveFR_run_cell),resolution_ratio);  
for i=1:1:length(aveFR_run_cell)
    aveFR_temp = aveFR_run_cell{i};
    mapping=linspace(1,length(aveFR_temp)+1,resolution_ratio);
    for j=1:1:length(aveFR_temp)
        aveFR_topo_mat(i,(mapping>=j)) = aveFR_temp(j);
    end 
end
aveFR_topo = mean(aveFR_topo_mat);
ste_aveFR_topo = std(aveFR_topo_mat,0,1)/sqrt(size(aveFR_topo_mat,1));

aveFR_topo_wholeWindow = [mean_aveFR_bfRun, aveFR_topo, mean_aveFR_aftRun];
steFR_topo_wholeWindow = [ste_aveFR_bfRun,ste_aveFR_topo, ste_aveFR_aftRun];

% plot
run_wholetime_aveFR_topo = figure;
subplot(2,1,1);
shadedErrorBar(x,aveFR_topo_wholeWindow,steFR_topo_wholeWindow);
ylabel ('firing rate');vline(time/30,'r'); vline((length(mean_speed_wholeWindow)-time)/30,'r');
xlim([0 10.33]); ylim([0 2]);
%xlim([0,length(mean_speed_wholeWindow)+5]);
%ylim([min(aveFR_topo_wholeWindow)-0.3,max(aveFR_topo_wholeWindow)+0.3]);
title(['FR topologized before, during, and after running deconvolve' sessions]);

subplot(2,1,2);
shadedErrorBar(x,speed_topo_wholeWindow*2*3.1415926*7.5/128,speed_ste_topo_wholeWindow*2*3.1415926*7.5/128);
ylabel ('speed (cm/s)');vline(time/30,'r'); vline((length(mean_speed_wholeWindow)-time)/30,'r');
xlim([0 10.33]); ylim([-5 20]);
xlabel('time(s)');
%xlim([0,length(mean_speed_wholeWindow)+5]);
%ylim([min(speed_topo_wholeWindow)-5,max(speed_topo_wholeWindow)+5]);
saveas(run_wholetime_aveFR_topo, [image_analysis_dest sessions, '_FR_run_allWindows_topo_deconvolve' num2str(threshold)]);

save([image_analysis_dest sessions '_spk_deconvolve_threshold' num2str(threshold) '.mat'],'aveFR_topo_wholeWindow','steFR_topo_wholeWindow','speed_topo_wholeWindow','speed_ste_topo_wholeWindow','-append');


%% topology - optimal resolution - test

% resolution_ratios  = [5,50,100,250,500,750,1000,1250,1500,1750,2000,2500,3000,4000,6000,8000,20000]; %divide the total process into how many parts
% 
% speed_topo_ave = zeros(1,length(resolution_ratios));
% for r = 1: length(resolution_ratios)
%     speed_topo_mat= zeros(length(speed_run_cell),resolution_ratios(r));
%     for i=1:1:length(speed_run_cell)
%         speed_temp = speed_run_cell{i};
%         mapping = linspace(1,length(speed_temp)+1,resolution_ratios(r));
%         for j = 1:1:length(speed_temp)
%             speed_topo_mat(i,(mapping>=j))=speed_temp(j);
%         end
%     end
%     speed_topo_ave(r) = std(mean(speed_topo_mat));
% end
% figure;plot(resolution_ratios,speed_topo_ave,'o'); % looks like numbers above 2000 is good
%scatter(resolution_ratios,speed_topo_ave);


