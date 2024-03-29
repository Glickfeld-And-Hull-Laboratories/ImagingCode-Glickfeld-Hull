%% SECTION - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
%sessions = {'190617_img1021_1','190617_img1023_1','190617_img1024_1',...
%    '190617_img1027_2','190618_img1025_1','200321_img1042_1','200321_img1049_1',...
%    '200321_img1063_1','200321_img1064_1'}; 
%days = {'1021-190617_1','1023-190617_1','1024-190617_1',...
%    '1027-190617_1','1025-190618_1','1042-200321_1','1049-200321_1',...
%    '1063-200321_1','1064-200321_1'};
sessions = {'190618_img1025_1'};
days = {'1025-190618_1'};
%there might be more than 1 sessions on a single subject on the same day
image_dest_base = 'Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\'; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'r','g','m','y','b'};

%% 1. when does the first peak happen?
for i = 1:length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{i}];
    behav_struct = load([behav_dest '\' days{i} '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    dfOvF_struct = load([image_dest, '_dfOvF_btmbaseline.mat']);
    dfOvF = dfOvF_struct.dfOvF_btmbase;
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    
    frames_runoff_mat = behav_struct.frames_runoff_mat;%frames*trials
    % generate a matrix for loooking at peak times after running offset
    offtime = 25; %2.5 seconds after running offset
    offset_mat = zeros(offtime,size(frames_runoff_mat,2)); %frames*trials
    for t = 1: size(frames_runoff_mat,2)
        offset_mat(:,t) = frames_runoff_mat(11,t):frames_runoff_mat(11,t)+offtime-1; % in frames_runoff_mat, the 11th frame in each trial is the first frame after running offset
    end
    dfOvF_offsetmat = zeros(size(offset_mat,2),size(offset_mat,1),size(dfOvF,1));%trials*frames*ROIs
    %x = -5:runTriggerDura-6;
    for n = 1:size(dfOvF,1)
        temp = dfOvF(n,:);
        dfOvF_offsetmat(:,:,n) = (temp(offset_mat))';
    end
    
    peak_loc = zeros(size(dfOvF_offsetmat,1),size(dfOvF,1));
    for t = 1:size(dfOvF_offsetmat,1) %for each trial
        for c = 1: size(dfOvF,1)
            [pks,locs] = findpeaks(dfOvF_offsetmat(t,:,c));
            %findpeaks(dfOvF_offsetmat(t,:,c)); %automatically plot the peaks
            %find all peaks --> if the first peak happens to be the biggest
            %one, then the first peak is the real peak.
            % if not, compare the peaks before the largest peak with the
            % largest one, if those peaks is not much smaller than the largest
            % one, then the previous peak is the first real peak.
            biggest = max(pks);
            inx = 1:length(locs);
            biggest_inx = inx(pks==biggest);
            if biggest_inx==1
                peak_loc(t,c) = locs(1);
            else
                for p = 1:biggest_inx
                    pk_diff = biggest - pks(p);
                    if pk_diff<0.05
                        peak_loc(t,c) = locs(p);
                        break; %once fullfills the requirement, jump out of the for loop
                    end
                end
            end
        end
    end
    
    %make a binary mask, only peak frame = 1, others = 0
    
    fig_dest = [image_dest_base sessions{i} '\runoffPeak\'];
    if ~exist(fig_dest)
        mkdir(fig_dest);
    end
    
    runoff_bi = zeros(size(dfOvF_offsetmat,1),offtime,size(dfOvF,1));%trial*frame*ROI
    for c = 1: size(dfOvF,1)
        for t = 1:size(dfOvF_offsetmat,1)
            runoff_bi(t,peak_loc(t,c),c)=1;
        end
        x = (0.1:0.1:2.5);
        y = [0 size(dfOvF_offsetmat,1)];
        figure; imagesc(x,y,runoff_bi(:,:,c));
        ylabel('trial #');
        xlabel('time(s)');
        title([sessions{i} ' ROI '  num2str(c)]);
        savefig([fig_dest sessions{i} '1peakTime_ROI' num2str(c)]);
    end
    save([image_dest '_imgAnalysis.mat'] ,'dfOvF_offsetmat','peak_loc','runoff_bi','-append'); 
end


%% 2. how many peaks does each trial have in total 
for i = 1:length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    dfOvF_offsetmat = img_anal.dfOvF_offsetmat;
    big_pks = zeros(size(dfOvF_offsetmat)); %trials*frames*ROIs
    sumpks = zeros(size(dfOvF_offsetmat,1),size(dfOvF_offsetmat,3));
    for t = 1:size(dfOvF_offsetmat,1) %for each trial
        for c = 1: size(dfOvF_offsetmat,3)
            [pks,locs] = findpeaks(dfOvF_offsetmat(t,:,c));
            ave_pks = mean(pks);
            for p = 1:length(pks)
                if pks(p)> ave_pks
                    big_pks(t,locs(p),c) = 1;
                end
            end
            sumpks(t,c) = sum(big_pks(t,:,c)==1);
        end
    end
    save([image_dest '_imgAnalysis.mat'] ,'big_pks','sumpks','-append');
end

%% is trial length and average speed correlated with when the first peak happens/how many peaks are there?
for i = 1:length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    peak_loc = img_anal.peak_loc;
    sumpks = img_anal.sumpks;
    behav_dest = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days{i}];
    behav_struct = load([behav_dest '\' days{i} '_behavAnalysis.mat']);
    speed = behav_struct.speed;
    frames_runoff_include = behav_struct.frames_runoff_include;
    trial_lens = cellfun(@length,frames_runoff_include);
    speed_runoff_include = cell(size(frames_runoff_include));
    for t = 1:size(frames_runoff_include,2)
        speed_runoff_include{t} = speed(frames_runoff_include{t});
    end
    avespd_runofftrials = cellfun(@mean,speed_runoff_include);
    
    fig_dest = [image_dest_base sessions{i} '\runoffPeak\'];
    %sumpks and ave speed/trial lengths
    for r = 1:size(sumpks,2)
        figure; scatter(avespd_runofftrials*2*3.1415926*7.5/128,sumpks(:,r),'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
        ylim([0 6]); %xlim([0 max(avespd_runofftrials*2*3.1415926*7.5/128)]);
        ylabel('total number of df/f peaks');
        xlabel('average speed (cm/s)');
        title([sessions{i} 'ROI' num2str(r)]);
        savefig([fig_dest sessions{i} 'sum_pks_Vs_speed_ROI' num2str(r)]);
        
        figure; scatter(trial_lens./10,sumpks(:,r),'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
        ylim([0 6]);
        ylabel('total number of df/f peaks');
        xlabel('running trial duration (s)');
        title([sessions{i} 'ROI' num2str(r)]);
        savefig([fig_dest sessions{i} 'sum_pks_Vs_duration_ROI' num2str(r)]);
    end
    
     for r = 1:size(peak_loc,2)
        figure; scatter(avespd_runofftrials*2*3.1415926*7.5/128,peak_loc(:,r),'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
        ylim([0 25]); %xlim([0 max(avespd_runofftrials*2*3.1415926*7.5/128)]);
        ylabel('time of first peak after running offset');
        xlabel('average speed (cm/s)');
        title([sessions{i} 'ROI' num2str(r)]);
        savefig([fig_dest sessions{i} '_1pk_time_Vs_speed_ROI' num2str(r)]);
        
        figure; scatter(trial_lens./10,peak_loc(:,r),'filled','MarkerFacecolor',[0.9373 0.396 0.2824]);
        ylim([0 25]);
        ylabel('time of first peak after running offset');
        xlabel('running trial duration (s)');
        title([sessions{i} 'ROI' num2str(r)]);
        savefig([fig_dest sessions{i} '_1pk_time_Vs_duration_ROI' num2str(r)]);
     end
    
end


