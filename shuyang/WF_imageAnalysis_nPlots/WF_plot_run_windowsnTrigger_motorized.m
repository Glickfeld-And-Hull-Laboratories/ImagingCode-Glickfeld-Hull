% plotting fluroscence data for wide field imaging, self-paced running experiments
% need: df/f and output from behav_analysis_movingDots_WF .
% generate cells for before and after running windows, 
% matrix for run trigger average and running windows, using findFrames_runWindows
% plot df/f and speed align to running onset and offset
% scatter plot for ave df/f before, during, and after running
% plot df/f from 300ms before running to 300ms after running, each point is the ave for 300ms

%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
% sessions = {'190617_img1021_1','190617_img1023_1','190617_img1024_1',...
%     '190618_img1025_1','200321_img1042_1','200321_img1049_1',...
%     '200321_img1064_1'}; 
% days = {'1021-190617_1','1023-190617_1','1024-190617_1',...
%     '1025-190618_1','1042-200321_1','1049-200321_1',...
%     '1064-200321_1'};

sessions = {'190617_img1021_1'};
days = {'1021-190617_1'};
image_dest_base    = 'Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\'; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'c','r','y','g'};


%% SECTION I: generate matrix of frames needed for plotting and load data needed for plotting
for ii = 1: length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_btmbaseline.mat']);
    dfOvF = dfOvF_strct.dfOvF_btmbase;
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    frames_run_cell = behav_struct.frames_run_cell;
    speed = behav_struct.speed;
    [frames_befo_run_cell,frames_aft_run_cell,frames_runTrigger_mat,...
    frames_runoff_include,frames_runoff_mat,frames_run_mat] = findFrames_runWindows (speed,frames_run_cell);
    
    save([behav_dest '\' days{ii} '_behavAnalysis.mat' ],...
        'frames_befo_run_cell','frames_aft_run_cell','frames_runTrigger_mat',...
        'frames_runoff_mat','frames_run_mat', 'frames_runoff_include','-append');
end


%%  SECTION IV: run trigger ave
for ii = 1: length(sessions)
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    speed = behav_struct.speed; speed = double(speed);
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_btmbaseline.mat']);
    dfOvF = dfOvF_strct.dfOvF_btmbase;
    frames_runTrigger_mat = behav_struct.frames_runTrigger_mat;

    runTrigSpeed = (speed(frames_runTrigger_mat))';
    ave_speed_runTrigger = mean(runTrigSpeed);
    ste_speed_runTrigger = std(runTrigSpeed,0,1)/sqrt(size(runTrigSpeed,1));
    
    % generate matrixes for df/f and plot----------------------------------
    dfOvF_runTrigger = zeros(size(frames_runTrigger_mat,2),size(frames_runTrigger_mat,1),size(dfOvF,1));
    for n = 1:size(dfOvF,1)
        temp = dfOvF(n,:);
        dfOvF_runTrigger(:,:,n) = (temp(frames_runTrigger_mat))'; %trials*frames*ROIs
%         set(0,'DefaultFigureVisible','off'); %doesn't show figures for the single sessions on the screen
%         windows_fig(n) = figure; 
%         %make every line start at the same height
%         startMean = mean(dfOvF_runTrigger(:,1,n));
%         diff = dfOvF_runTrigger(:,1,n)-startMean;
%         diff_array = repmat(diff,1, size(dfOvF_runTrigger,2));
%         dfOvF_runTrigger_plot = dfOvF_runTrigger(:,:,n) - diff_array;
%         dfOvF_runTrigger_plot = dfOvF_runTrigger_plot';
%         plot(x,dfOvF_runTrigger_plot);
%         xlabel('frames');
%         ylabel('df/f');
%         set(gca,'XTick',x);
%         vline(-1, 'k','running start');
%         title(['df/f for each running window ROI' num2str(n)]); 
%         saveas(windows_fig(n), [image_dest '_runTrigSessions_ROI' num2str(n)]);
    end
    
    ave_dfOvF_runTrigger = squeeze(mean(dfOvF_runTrigger,1)); %average across trials
    ste_dfOvF_runTrigger = squeeze(std(dfOvF_runTrigger,0,1)/sqrt(size(dfOvF_runTrigger,1)));
    
    ave_dfOvF_runTrigger = ave_dfOvF_runTrigger';
    ste_dfOvF_runTrigger = ste_dfOvF_runTrigger';
    
    set(0,'DefaultFigureVisible','on'); 
    x = ((0:length(ave_dfOvF_runTrigger)-1)/10)-1;
    if size(ave_dfOvF_runTrigger,2)>1
        xplot = repmat(x,size(ave_dfOvF_runTrigger,1),1);
    else
        xplot = x;
    end
    colorcode = {[0.4000 0.7608 0.6471],[0.9882 0.5529 0.3843],[0.5529 0.6275 0.7961],[0.9059 0.5412 0.7647]};%,[0.6510 0.8471 0.3294]}; %brewermap set2
    eg_runonset = figure; 
    %shadedErrorBar(x,dfOvF_runtrig_slow_session,ste_runtrig_slow_cells,{'color',[0.1373
    %0.5451 0.2706]}); it's ugly use shadederrorbar and need a for loop, shadederrobar only takes vector
    subplot(2,1,1);hold on;
    for r = 1:size(ave_dfOvF_runTrigger,1)
        errorbar(xplot(r,:),ave_dfOvF_runTrigger(r,:),ste_dfOvF_runTrigger(r,:),...
            '.','LineStyle','-','linewidth', 1,'MarkerSize',4,'Color',colorcode{r});hold on;
    end
    xlim([-1.1 1.6]);
    ylim([-0.1 0.45]);
    vline(0,'k');
    %legend('ROI1','ROI2','ROI3','ROI4');legend('boxoff'); 
    ylabel('df/F');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',7);
    
    subplot(2,1,2);hold on;
    errorbar(x,ave_speed_runTrigger*2*3.1415926*7.5/128,ste_speed_runTrigger*2*3.1415926*7.5/128,...
        '.','LineStyle','-','linewidth', 1,'MarkerSize',4,'Color',[0.4500 0.4500 0.4500]);hold on;
    xlim([-1.1 1.6]);
    xlabel('time from running onset(s)');
    ylabel('speed(cm/s)');
    vline(0, 'k');
    ylim([0 20]);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',7);
    %supertitle(['run triggered average',sessions{ii}] );
    %saveas(mean_fig, [image_dest '_runTrigAve']);
    eg_runonset.Units = 'centimeters';
    eg_runonset.Position = [3 3 5.5 5];
    fig_name = ['eg_session_runonset_',sessions{ii}];
    path = 'Z:\Analysis\figures\figure1_WF\';
    orient(eg_runonset,'landscape')
    print(eg_runonset,[path,fig_name],'-r600','-depsc');
    %saveas(eg_runonset, [image_dest '_runoffAve']);
    %save([image_dest '_imgAnalysis.mat' ],'ave_dfOvF_runTrigger', 'ave_speed_runTrigger', '-append');
    
end


%%  SECTION V: run offset ave
for ii = 1: length(sessions)
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    speed = behav_struct.speed; speed = double(speed);
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_btmbaseline.mat']);
    dfOvF = dfOvF_strct.dfOvF_btmbase;
    frames_runoff_mat = behav_struct.frames_runoff_mat;

    runOffSpeed = (speed(frames_runoff_mat))';
    ave_speed_runoff = mean(runOffSpeed);
    ste_speed_runoff = std(runOffSpeed,0,1)/sqrt(size(runOffSpeed,1));
    
    % generate matrixes for df/f and plot----------------------------------
    dfOvF_runOff = zeros(size(frames_runoff_mat,2),size(frames_runoff_mat,1),size(dfOvF,1));%trials*frames*ROIs
    %x = -5:runTriggerDura-6;
    for n = 1:size(dfOvF,1)
        temp = dfOvF(n,:);
        dfOvF_runOff(:,:,n) = (temp(frames_runoff_mat))';
%         set(0,'DefaultFigureVisible','off'); %doesn't show figures for the single sessions on the screen
%         windows_fig(n) = figure; 
%         %make every line start at the same height
%         startMean = mean(dfOvF_runOff(:,1,n));
%         diff = dfOvF_runOff(:,1,n)-startMean;
%         diff_array = repmat(diff,1, size(dfOvF_runOff,2));
%         dfOvF_runOff_plot = dfOvF_runOff(:,:,n) - diff_array;
%         plot(x,dfOvF_runOff_plot');
%         xlabel('frames');
%         ylabel('df/f');
%         set(gca,'XTick',x);
%         vline(-1, 'k','running start');
%         title(['df/f for each running window ROI' num2str(n)]); 
%         saveas(windows_fig(n), [image_dest '_runTrigSessions_ROI' num2str(n)]);
    end
    
    ave_dfOvF_runOff = squeeze(mean(dfOvF_runOff,1)); %average across trials
    ste_dfOvF_runOff = squeeze(std(dfOvF_runOff,0,1)/sqrt(size(dfOvF_runOff,1)));
    
    ave_dfOvF_runOff = ave_dfOvF_runOff';
    ste_dfOvF_runOff = ste_dfOvF_runOff';
    
    set(0,'DefaultFigureVisible','on');
    x = ((0:length(ave_dfOvF_runOff)-1)/10)-1;
    if size(ave_dfOvF_runOff,2)>1
        xplot = repmat(x,size(ave_dfOvF_runOff,1),1);
    else
        xplot = x;
    end
    
    x = ((0:length(ave_dfOvF_runOff)-1)/10)-1;
    mean_fig = figure; 
    subplot(2,1,1);hold on;
    for r = 1:size(ave_dfOvF_runOff,1)
        errorbar(xplot(r,:),ave_dfOvF_runOff(r,:),ste_dfOvF_runOff(r,:),...
            '.','LineStyle','-','linewidth', 1,'MarkerSize',4,'Color',colorcode{r});hold on;
    end
    xlim([-1.1 1.6]);
    ylim([-0.1 0.45]);
    %legend({'ROI1','ROI2','ROI3','ROI4'},'FontSize',10,'location','southeast');legend('boxoff');
    vline(0,'k');
    ylabel('df/F');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',7);
    %ax = gca;
    %ax.LabelFontSizeMultiplier = 1.25;
    
    subplot(2,1,2);hold on;
    errorbar(x,ave_speed_runoff*2*3.1415926*7.5/128,ste_speed_runoff*2*3.1415926*7.5/128,...
        '.','LineStyle','-','linewidth', 1,'MarkerSize',4,'Color',[0.45 0.45 0.45]);hold on;
    xlabel('time from running offset(s)');
    ylabel('speed(cm/s)');
    xlim([-1.1 1.6]);
    ylim([0 20]);
    vline(0, 'k');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',7);
    %ax = gca;
    %ax.LabelFontSizeMultiplier = 1.25;
    
    %supertitle(['run triggered average',sessions{ii}] );
    %saveas(mean_fig, [image_dest '_runTrigAve']);
    mean_fig.Units = 'centimeters';
    mean_fig.Position = [1 1 5.5 5];
    fig_name = ['eg_session_runoffset_',sessions{ii}];
    path = 'Z:\Analysis\figures\figure1_WF\';
    orient(mean_fig,'landscape');
    print(mean_fig,[path,fig_name],'-r600','-depsc');
   
    %supertitle(['run offset average',sessions{ii}] );
%     saveas(mean_fig, [image_dest '_runoffAve']);
%     save([image_dest '_imgAnalysis.mat' ],'ave_dfOvF_runOff', 'ave_speed_runoff', '-append');

    
end


%% SECTION II:draw scatter plot for df/f of right before and right after every running window
for ii = 1:length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_btmbaseline.mat']);
    dfOvF = dfOvF_strct.dfOvF_btmbase;
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    frames_run_cell = behav_struct.frames_run_cell; 
    frames_run = cell2mat(frames_run_cell);
    frames_befo_run = behav_struct.frames_befo_run_cell;
    frames_befo_run = cell2mat(frames_befo_run);
    frames_aft_run = behav_struct.frames_aft_run_cell;
    frames_aft_run = cell2mat(frames_aft_run);
    speed = behav_struct.speed; speed = double(speed);

    
    % calculate ave of dfOvF and put before, during run, after run together-----------------------------------------------
    dfOvF_befoRun = dfOvF(:,frames_befo_run);
    % this give you dfOvF of first, second, and third column in frames_befor_run.
    %ROI1 are now placed in a row(1-78:first column of frames_befo_run; 79-157: second column; 158-234: third column)
    dfOvF_aftRun = dfOvF(:,frames_aft_run);
    dfOvF_run = dfOvF(:,frames_run);
    
    ave_dfOvF_befoRun = mean(dfOvF_befoRun');
    %tranpose it to each column is an ROI later
    ste_dfOvF_befoRun = std(dfOvF_befoRun')/sqrt(length(dfOvF_befoRun'));
    ave_dfOvF_run = mean(dfOvF_run');
    ste_dfOvF_run = std(dfOvF_run')/sqrt(length(dfOvF_run'));
    ave_dfOvF_aftRun = mean(dfOvF_aftRun');
    ste_dfOvF_aftRun = std(dfOvF_aftRun')/sqrt(length(dfOvF_aftRun'));
    ave_befoRunaft = [ave_dfOvF_befoRun',ave_dfOvF_run',ave_dfOvF_aftRun'];
    ste_befoRunaft = [ste_dfOvF_befoRun',ste_dfOvF_run',ste_dfOvF_aftRun'];
    
    %do the same thing for speed-------------------------------------------
    speed_befoRun = speed(frames_befo_run);
    speed_aftRun = speed(frames_aft_run);
   
    ave_speed_befoRun = mean(speed_befoRun);
    %tranpose it to each column is an ROI later
    ste_speed_befoRun = std(speed_befoRun)/sqrt(length(speed_befoRun));
    ave_speed_aftRun = mean(speed_aftRun);
    ste_speed_aftRun = std(speed_aftRun)/sqrt(length(speed_aftRun));
    
    %----------------------------------------------------------------------
    dfOvF_run_buffer = figure;
    x = [1,2,3]; x_plot = repmat(x,size(dfOvF,1),1);
    errorbar(x_plot',ave_befoRunaft',ste_befoRunaft','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); 
    legend('ROI1','ROI2','ROI3','ROI4');
    xlim([0.5 3.5]);
    %xlabel ('');
    x1= [1,2,3];
    set(gca,'XTick',x1,'XTicklabel',{'right before','run','right after'});
    ylabel('df/f');
    title(['df/f right before and after running',sessions{ii}]); axis square;
    saveas(dfOvF_run_buffer, [image_dest '_dfOvF_runVsSurround']);
    save([behav_dest '\' days{ii} '_behavAnalysis.mat' ],...
        'ave_dfOvF_befoRun','ave_dfOvF_aftRun','ste_dfOvF_befoRun','ste_dfOvF_aftRun',...
        'ave_speed_befoRun','ste_speed_befoRun','ave_speed_aftRun','ste_speed_aftRun','-append');
    save([image_dest '_imgAnalysis.mat' ],'ave_befoRunaft','-append');
end

%{
%% SECTION III: df/f before, after and during running, every 300ms.
for ii = 1: length(sessions)
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{ii}];
    behav_struct = load([behav_dest '\' days{ii} '_behavAnalysis.mat']);
    speed = behav_struct.speed; speed = double(speed);
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    dfOvF_strct = load([image_dest, '_dfOvF_btmbaseline.mat']);
    dfOvF = dfOvF_strct.dfOvF_btmbase;
    frames_run_mat = behav_struct.frames_run_mat;
    
    %frames_run_buffer_mat = behav_output.frames_run_buffer_mat;
    ave_dfOvF_befoRun = behav_struct.ave_dfOvF_befoRun;
    ave_dfOvF_aftRun = behav_struct.ave_dfOvF_aftRun;
    ste_dfOvF_befoRun = behav_struct.ste_dfOvF_befoRun;
    ste_dfOvF_aftRun = behav_struct.ste_dfOvF_aftRun;
    ave_speed_befoRun = behav_struct.ave_speed_befoRun;
    ste_speed_befoRun = behav_struct.ste_speed_befoRun;
    ave_speed_aftRun = behav_struct.ave_speed_aftRun;
    ste_speed_aftRun = behav_struct.ste_speed_aftRun;
    
    % calculate average speed of every 300ms ------------------------------
    speed_run_mat = nan(size(frames_run_mat,1), size(frames_run_mat,2));
    for v = 1:size(speed_run_mat,1)
        temp = frames_run_mat(v,~isnan(frames_run_mat(v,:)));
        speed_run_mat(v,1:size(temp,2)) = speed(temp);
    end
    speed_run_mat(speed_run_mat(v)==0)=nan;
    mean_columns = mean(speed_run_mat,'omitnan');
    bin = floor(length(mean_columns)/3);
    mean_run_every300ms = zeros(1,bin);
    ste_run_every300ms = zeros(1,bin);
    for x = 1:bin
        y = 3*x -2;
        mean_run_every300ms(x) = mean(mean_columns(y:y+2));
        ste_run_every300ms(x) = std(mean_columns(y:y+2))/sqrt(3);
    end
   % if length(mean_columns) can't be divided by 3 in total, you're throwing away the remainers.
   
    % calculate average fluorescence of every 300ms------------------------
    dfOvF_staybase_mat = nan(size(frames_run_mat,1), size(frames_run_mat,2));
    dfOvF_run_cell = {};
    mean_dfOvF_300ms = []; ste_dfOvF_300ms = [];
    for q = 1:size(dfOvF,1)
        for n = 1: size(dfOvF_staybase_mat,1)
            temp = frames_run_mat(n,~isnan(frames_run_mat(n,:)));
            dfOvF_run_cell{q}(n,1:size(temp,2)) = dfOvF(q,temp);
        end
        dfOvF_run_cell{q}(dfOvF_run_cell{q}==0)=nan;
        mean_columns = mean(dfOvF_run_cell{q},'omitnan');
        %calculate average for every 300ms
        bin = floor(length(mean_columns)/3);
        mean_every300ms = zeros(1,bin);
        ste_every300ms = zeros(1,bin);
        for x = 1:bin
            y = 3*x -2;
            mean_every300ms(x) = mean(mean_columns(y:y+2));
            ste_every300ms(x) = std(mean_columns(y:y+2))/sqrt(3);
        end
       mean_dfOvF_300ms(q,:) = mean_every300ms;
       ste_dfOvF_300ms(q,:) = ste_every300ms;
    end
        
    % plot ROIs and speed--------------------------------------------------
    %plotting the 1st and last 900ms of every running window, running
    %period in between is omitted!!!!!!!!!!!!!!!!!!!
    dfOvF_plot300ms = [ave_dfOvF_befoRun',mean_dfOvF_300ms,ave_dfOvF_aftRun'];
    dfOvF_plot300ms = [dfOvF_plot300ms(:,1:10),dfOvF_plot300ms(:,end-9:end)];
    ste_dfOvF_plot300ms = [ste_dfOvF_befoRun',ste_dfOvF_300ms,ste_dfOvF_aftRun'];
    ste_dfOvF_plot300ms = [ste_dfOvF_plot300ms(:,1:10),ste_dfOvF_plot300ms(:,end-9:end)];
    speed_plot300ms = [ave_speed_befoRun,mean_run_every300ms,ave_speed_aftRun];
    speed_plot300ms = [speed_plot300ms(1:10),speed_plot300ms(end-9:end)];
    ste_speed_plot300ms = [ste_speed_befoRun,ste_run_every300ms,ste_speed_aftRun];
    ste_speed_plot300ms = [ste_speed_plot300ms(1:10),ste_speed_plot300ms(end-9:end)];
    
    dfOvF_run_300ms = figure;
    subplot(2,1,1);
    errorbar(dfOvF_plot300ms',ste_dfOvF_plot300ms','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20);legend
    ylabel ('df/f');

    subplot(2,1,2);
    errorbar(speed_plot300ms, ste_speed_plot300ms,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20);
    ylabel('speed');
    %xticks([ ])
    %xticklabels({'','','','','',''});
    supertitle(['df/f before, during, and after running',sessions{ii}]);
    saveas(dfOvF_run_300ms, [image_dest, '_dfOvF_runWindows']);
     
    %set(gca,'XTick',x,'XTicklabel',{'before run','1','','','','5','','','','',...
        %'10','','','','','15','','','','','20','','after run'});
    %'','','','25','','','', '','30',...
    % '','','','','35','','','','','40','','','','','45','','','','','50','','',...
    %   '','','55','','','','','60','','','','','65','','','','','70','','','','',...
    %  '75','','','','','80','','','','','85','','','','','90','','','','','','',...
    %  '','','','','','','','','','','','','','','','','',
 
    % legend(handle,'ROI1', 'ROI2');
    save([image_dest '_imgAnalysis.mat' ],'dfOvF_plot300ms', 'speed_plot300ms', '-append');

end
%}

