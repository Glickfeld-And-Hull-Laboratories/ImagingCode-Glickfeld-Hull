% does the amplitude of Ca transient change under different conditions?
% condition 1 : stationary (used to be not including 1s after running offset and 1s before running onset, but now includes whole window)
% condition 2: running (used to be not including 1s after running onset and 1s? before running offset, but now includes whole window)
% condition 3: airpuff response during running (within 200ms of airpuff onset)
% condition 4: airpuff response during stationary (within 200ms of airpuff onset)
% events need to be isolated: no other events within 600ms before and/or after
% frm_stay_midpart: all stationary frames without beginning and end of stationary
% isolated_inx_neurons: index of isolated event peaks 1*#of cells, each cell element is a vector
% spk_iso_stay_neurons: index of isolated event peaks during stationary, 1*# of neurons, each cell element is a vector
% isoevent_stay_neurons: index of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
% dfOvF_btm_stay_isoevents: %df/f of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
% dfOvF_btm_stay_isoevents_neurons: average df/f of isolated event of each cell, averaged across events, neurons*frames
% ave_dfOvF_btm_stay_iso_session: average df/f of isolated event across cells, use this to plot

%% Section I: set paths and create analysis folders for each session
%define the directory and files
clear;
sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};
%sessionID = {'1023-190923_1','','','','','','','',''};
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; 
color_code = {'c','r','y','g'};

%% SectionII: for each session: running experiment, event amplitude during running vs stationary
% get events from each session, and save all of them into one variable at
% the end of the for loop
for ii = 1: length(sessions)
   % load data
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\']; 
    rawF_output = load([image_analysis_dest sessions{ii},'_deconvolution_thresh-4_TCave_cl.mat']);
    rawF = rawF_output.TCave_cl;
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{ii} '\'];
    behav_output = load([behav_dest days{ii} '_behavAnalysis.mat']);
    %behav_output = load([behav_dest days{i} '_first18000frames_behavAnalysis.mat']);
    frm_stay_cell = behav_output.frames_stay_cell;
    frm_run_cell = behav_output.frames_run_cell;
    frm_stay = cell2mat(frm_stay_cell);
    frm_run = cell2mat(frm_run_cell);
    frm_run_cell = behav_output.frames_run_cell;
    frm_run_starts = zeros(1,size(frm_run_cell,2));
    for r = 1:size(frm_run_cell,2)
        frm_run_starts(r) = frm_run_cell{r}(1);
    end
    % stationary: not including 1s (30 frames) after running offset and 1s before
    frm_stay_midpart_cell = {};
    for t = 1: size(frm_stay_cell,2)
        if length(frm_stay_cell{t})>= 61
        frm_stay_midpart_cell = [frm_stay_midpart_cell frm_stay_cell{t}(31:end-30)];
        end
    end
    frm_stay_midpart = cell2mat(frm_stay_midpart_cell); % all stationary frames without beginning and end of stationary
    
    %running (not including 1s after running onset and 1s? before running offset)
    frm_run_midpart_cell = {};
    for t = 1: size(frm_run_cell,2)
        if length(frm_run_cell{t})>= 61
        frm_run_midpart_cell = [frm_run_midpart_cell frm_run_cell{t}(31:end-30)];
        end
    end
    if isempty(frm_run_midpart_cell) == 0
        frm_run_midpart = cell2mat(frm_run_midpart_cell);
    end
    
    % for each cell,find isolated events in running and stationary (no event within 650ms before that event)
    threshold = -4;
    spk_deconv_output = load([image_analysis_dest sessions{ii},'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
    spk_inx_neurons = spk_deconv_output.spk_inx_cl;
    
    % identify events that are isolated: 600ms = 18 frames. need to do this then see what behavioral state they're in, not the other way because
    % if you find intersect of peak index and behavioral state frame first, then the diff between the frames you find can be big if there's a different behavioral state in between.
    noother_after_cells = cell(1,size(spk_inx_neurons,2));
    noother_befo_cells = cell(1,size(spk_inx_neurons,2));
    isolated_inx_neurons = cell(1,size(spk_inx_neurons,2));
    for c = 1: size(spk_inx_neurons,2)
        noother_after_inx = find(diff(spk_inx_neurons{c})>18);% this means the indexed frame# being found here has a difference of >18 with the next frame#
        noother_after_cells{c} = spk_inx_neurons{c}(noother_after_inx);
        noother_befo_inx = noother_after_inx+1; % the diff of inx-(inx-1) > 18
        noother_befo_cells{c} = spk_inx_neurons{c}(noother_befo_inx);
        isolated_inx_neurons{c} = intersect(noother_after_cells{c},noother_befo_cells{c}); %index of isolated event peaks 1*#of cells, each cell element is a vector
    end
    
    %isolated event peaks and matrices during stationary and running
    spk_iso_stay_neurons = cell(1,size(spk_inx_neurons,2));% index of isolated event peaks during stationary, 1*# of neurons, each cell element is a vector
    spk_iso_run_neurons = cell(1,size(spk_inx_neurons,2));
    isoevent_stay_neurons = cell(1, size(spk_inx_neurons,2));% index of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
    isoevent_run_neurons = cell(1, size(spk_inx_neurons,2));
    rawF_stay_isoevents = cell(1, size(spk_inx_neurons,2)); %df/f of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
    rawF_run_isoevents = cell(1, size(spk_inx_neurons,2));
    rawF_stay_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);% average df/f of isolated event of each cell, averaged across events, neurons*frames
    rawF_run_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);
    rawF_spk_run_begin_neurons = zeros(size(spk_inx_neurons,2),30);
    rawF_spk_run_middle_neurons = zeros(size(spk_inx_neurons,2),30);
    rawF_spk_run_late_neurons = zeros(size(spk_inx_neurons,2),30);
    rawF_stay_rand_isoevents = [];
    rawF_run_rand_isoevents = [];
    nevent_stay = [];
    nevent_run = [];
    nevent = zeros(1,size(spk_inx_neurons,2));
    neventrun = zeros(1,size(spk_inx_neurons,2));
    spk_iso_run_neurons_vec = [];
    isospk_run_begin = 0;
    isospk_run_middle = 0;
    isospk_run_late = 0;
    rawF_spk_run_begin = cell(1,size(spk_inx_neurons,2));
    rawF_spk_run_middle = cell(1,size(spk_inx_neurons,2));
    rawF_spk_run_late = cell(1,size(spk_inx_neurons,2));
    for c = 1: size(spk_inx_neurons,2)% for each cell
        %stationary
        spk_iso_stay_neurons{c} = intersect(frm_stay,isolated_inx_neurons{c});
        %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
        isoevent_stay_neurons{c} = zeros(length(spk_iso_stay_neurons{c}),30); % for each cell, create a matrix, # of events*30(1s), events*frames, each row is a stationary event
        rawF_stay_isoevents{c} =  zeros(length(spk_iso_stay_neurons{c}),30);% number of events*frames
        rawF_c = rawF(:,c);
        for s = 1:length(spk_iso_stay_neurons{c})
             isoevent_stay_neurons{c}(s,:) =  spk_iso_stay_neurons{c}(s)-11:spk_iso_stay_neurons{c}(s)+18;
             rawF_stay_isoevents{c}(s,:) = rawF_c(isoevent_stay_neurons{c}(s,:));
        end
        
        %running
        if isempty(frm_run_midpart_cell) == 0 % if there are running trials that are longer than 2s, then continue
            spk_iso_run_neurons{c} = intersect(frm_run,isolated_inx_neurons{c});
            spk_iso_run_neurons_vec = cat(1, spk_iso_run_neurons_vec,spk_iso_run_neurons{c});
            %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
            isoevent_run_neurons{c} = zeros(length(spk_iso_run_neurons{c}),30);
            rawF_run_isoevents{c} =  zeros(length(spk_iso_run_neurons{c}),30); %each row is a spike, each column is a frame
            for s = 1:length(spk_iso_run_neurons{c}) % for each spike of each cell
                isoevent_run_neurons{c}(s,:) =  spk_iso_run_neurons{c}(s)-11:spk_iso_run_neurons{c}(s)+18;
                rawF_run_isoevents{c}(s,:) = rawF_c(isoevent_run_neurons{c}(s,:));
                %-------------------------------------------------------------------------------
                % decide when the spike happens during running and
                % catagorize them based on when it happens --- for comparing the amplitude of isolated events during different running period
                % first find the first frame of that running window,
                a = frm_run_starts(frm_run_starts <= spk_iso_run_neurons{c}(s)); % these are incides of all first frames of running that is before the spike
                tfromstart = spk_iso_run_neurons{c}(s) - a(end); %the last one in a is the first frame that is in the same running window as the spike
                if tfromstart <= 15
                    isospk_run_begin = isospk_run_begin + 1;
                    rawF_spk_run_begin{c} = cat(1,rawF_spk_run_begin{c},rawF_run_isoevents{c}(s,:));
                elseif tfromstart > 15 && tfromstart <= 30
                    isospk_run_middle = isospk_run_middle + 1;
                    rawF_spk_run_middle{c} = cat(1,rawF_spk_run_middle{c},rawF_run_isoevents{c}(s,:));
                else
                    isospk_run_late = isospk_run_late + 1;
                    rawF_spk_run_late{c} = cat(1,rawF_spk_run_late{c},rawF_run_isoevents{c}(s,:));
                end
            end
            
            % draw same number of events from stationary and running and then do average-----------------------------------------------------
            if isempty(rawF_run_isoevents{c}) == 0 % for each cell, if there is isolated spike during running
                nevent_stay = size(rawF_stay_isoevents{c},1);
                nevent_run = size(rawF_run_isoevents{c},1);
                if nevent_stay >= nevent_run % if there're more events during stationary than running
                    nevent(c) = nevent_run;
                    %randomly draw same number of events during stationary as during running
                    rand_stay = randperm(nevent_stay,nevent(c));
                    rawF_stay_rand_isoevents = rawF_stay_isoevents{c}(rand_stay,:);
                    % average across trials
                    rawF_stay_isoevents_neurons(c,:) = mean(rawF_stay_rand_isoevents);
                    % average across trials
                    rawF_run_isoevents_neurons(c,:) = mean(rawF_run_isoevents{c});
                else % if there're more events during running
                    nevent(c) = nevent_stay;
                    rand_run = randperm(nevent_run,nevent(c));
                    rawF_run_rand_isoevents = rawF_run_isoevents{c}(rand_run,:);
                    rawF_stay_isoevents_neurons(c,:) = mean(rawF_stay_isoevents{c});
                    rawF_run_isoevents_neurons(c,:) = mean(rawF_run_rand_isoevents);
                end
            end
            
            %for running, draw same number of events from different running times and then do average
            if isempty(rawF_spk_run_begin{c}) == 0 && isempty(rawF_spk_run_middle{c}) == 0 && isempty(rawF_spk_run_late{c}) == 0 % for each cell, if there is isolated spike during every running period
                nevent_begin = size(rawF_spk_run_begin{c},1);
                nevent_middle = size(rawF_spk_run_middle{c},1);
                nevent_late = size(rawF_spk_run_late{c},1);
                neventrun(c) = min([nevent_begin,nevent_middle,nevent_late]);
                rand_begin = randperm(nevent_begin,neventrun(c));
                rand_middle = randperm(nevent_middle,neventrun(c));
                rand_late = randperm(nevent_late,neventrun(c));
                rawF_spk_run_begin_rand = rawF_spk_run_begin{c}(rand_begin,:);
                rawF_spk_run_begin_neurons(c,:) = mean(rawF_spk_run_begin_rand); % average across trials
                rawF_spk_run_middle_rand = rawF_spk_run_middle{c}(rand_middle,:);
                rawF_spk_run_middle_neurons(c,:) = mean(rawF_spk_run_middle_rand);
                rawF_spk_run_late_rand = rawF_spk_run_late{c}(rand_late,:);
                rawF_spk_run_late_neurons(c,:) = mean(rawF_spk_run_late_rand);
            end
        end
    end
    
    % if there's a cell that doesn't spike during running,there will be a line of NaNs in dfOvF_btm_run/stay_isoevents_neurons, need to delete this Nan line before doing average
    rawF_run_isoevents_neurons = rawF_run_isoevents_neurons(all(~isnan(rawF_run_isoevents_neurons),2),:);
    ave_rawF_run_iso_session = mean(rawF_run_isoevents_neurons);% average across cells
    ste_rawF_run_iso_session = std(rawF_run_isoevents_neurons,0,1)/sqrt(size(rawF_stay_isoevents_neurons,1)); % ste
    
    rawF_stay_isoevents_neurons = rawF_stay_isoevents_neurons(all(~isnan(rawF_stay_isoevents_neurons),2),:);
    ave_rawF_stay_iso_session = mean(rawF_stay_isoevents_neurons);
    ste_rawF_stay_iso_session = std(rawF_stay_isoevents_neurons,0,1)/sqrt(size(rawF_stay_isoevents_neurons,1));
    
    % for dfOvF_spk_run_begin/middle/late_neurons, delete those lines with zeros, that means that cell doesn't fire during all 3 conditions
    rawF_spk_run_begin_neurons = rawF_spk_run_begin_neurons(any(rawF_spk_run_begin_neurons~=0,2),:);
    rawF_spk_run_middle_neurons = rawF_spk_run_middle_neurons(any(rawF_spk_run_middle_neurons~=0,2),:);
    rawF_spk_run_late_neurons = rawF_spk_run_late_neurons(any(rawF_spk_run_late_neurons~=0,2),:);
    
%     %count how many events you have in total
%     nevents_run = sum(cellfun(@(x)x(1),cellfun(@size,isoevent_run_neurons,'UniformOutput',false)));% # of isolated events for all cells during running
%     nevents_stay = sum(cellfun(@(x)x(1),cellfun(@size,isoevent_stay_neurons,'UniformOutput',false)));
%     
    %plot
    x = (1:30)/30;
    Caevent = figure;
    errorbar(x,ave_rawF_stay_iso_session,ste_rawF_stay_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
    errorbar(x,ave_rawF_run_iso_session,ste_rawF_run_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
    % number of calcium events is wrong in the current legend
    legend('stationary', 'running');
    ylabel('rawF');
    xlabel('time(s)');
    title(sessions(ii));
    text(0.75,0.45, ['nevents = ' num2str(sum(nevent))]);
    text(0.75,0.4, ['nPCs = ' num2str(size(rawF_stay_isoevents_neurons,1))]);
    saveas(Caevent,[image_analysis_dest '\' days{ii} '_CaAmp']);
    % save variables
    
    Caevent_run = figure;
    fast_errbar(x,rawF_spk_run_begin_neurons,1,'color',[0.1373 0.5451 0.2706]);hold on;
    fast_errbar(x,rawF_spk_run_middle_neurons,1,'color',[0.2549 0.6706 0.3647]); hold on;
    fast_errbar(x,rawF_spk_run_late_neurons,1,'color',[0.7294 0.8941 0.7020]); hold on;
    legend('0-0.5s','0.5-1s','later than 1s');
    ylabel('rawF');
    xlabel('time(s)');
    title(sessions(ii));
    text(0.75,0.4, ['nPCs = ' num2str(size(rawF_spk_run_begin_neurons,1))]);
    saveas(Caevent_run,[image_analysis_dest '\' days{ii} '_CaAmp_run']);
    
    % if isempty(frm_run_midpart_cell) == 0
    save([image_analysis_dest sessions{ii} '_isoCaEvent.mat'],'rawF_stay_isoevents','rawF_run_isoevents',...
        'rawF_stay_isoevents_neurons','rawF_run_isoevents_neurons',...
        'ave_rawF_stay_iso_session','ste_rawF_stay_iso_session',...
        'ave_rawF_run_iso_session','ste_rawF_run_iso_session','nevent',...
        'rawF_spk_run_begin_neurons','rawF_spk_run_middle_neurons','rawF_spk_run_late_neurons',...
        '-append');
    % else
%     fprintf(['no running trials longer than 2s in this session ' sessions{ii}]);
%     save([image_analysis_dest sessions{ii} '_isoCaEvent.mat'],'frm_stay_midpart',...
%         'frm_run_midpart','isolated_inx_neurons','spk_iso_stay_neurons',...
%         'isoevent_stay_neurons','dfOvF_btm_stay_isoevents','dfOvF_btm_stay_isoevents_neurons',...
%         'ave_dfOvF_btm_stay_iso_session','ste_dfOvF_btm_stay_iso_session','nevents_stay');
end
    


