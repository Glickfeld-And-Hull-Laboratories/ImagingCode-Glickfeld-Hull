% does the amplitude of Ca transient change under different conditions?
% condition 1 : stationary 
% condition 2: fast running 
% condition 3: slow running
% events need to be isolated: no other events within 600ms before and/or after
% isolated_inx_neurons: index of isolated event peaks 1*#of cells, each cell element is a vector
% spk_iso_stay_neurons: index of isolated event peaks during stationary, 1*# of neurons, each cell element is a vector
% isoevent_stay_neurons: index of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
% dfOvF_btm_stay_isoevents: %df/f of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
% dfOvF_btm_stay_isoevents_neurons: average df/f of isolated event of each cell, averaged across events, neurons*frames
% ave_dfOvF_btm_stay_iso_session: average df/f of isolated event across cells, use this to plot

%% Section I: set paths and create analysis folders for each session
%define the directory and files

clear;
sessions = {'200116_img1041','200214_img1042','200217_img1061','200225_img1049','200319_img1064','200319_img1064_2'};
days = {'1041-200116_1','1042-200214_1','1061-200217_1','1049-200225_1','1064-200319_1','1064-200319_2'};
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\'; 
color_code = {'c','r','y','g'};

%% SectionII: for each session: running experiment, event amplitude during running vs stationary
% get events from each session, and save all of them into one variable at
% the end of the for loop
for ii = 1: length(sessions)
    % load data
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvF_btm = dfOvF_strct.dfOvF_btm_cl;
    behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\running\behavioral_analysis\' days{ii} '\'];
    behav_output = load([behav_dest days{ii} '_behavAnalysis.mat']);
    run_fast_vec = behav_output.run_fast_vec;
    run_slow_vec = behav_output.run_slow_vec;
    stay_vec = behav_output.stay_vec;
    
    % for each cell,find isolated events in running and stationary (no event within 650ms before that event)
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
    spk_iso_runfast_neurons = cell(1,size(spk_inx_neurons,2));
    spk_iso_runslow_neurons = cell(1,size(spk_inx_neurons,2));
    isoevent_stay_neurons = cell(1, size(spk_inx_neurons,2));% index of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
    isoevent_runfast_neurons = cell(1, size(spk_inx_neurons,2));
    isoevent_runslow_neurons = cell(1, size(spk_inx_neurons,2));
    dfOvF_btm_stay_isoevents = cell(1, size(spk_inx_neurons,2)); %df/f of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
    dfOvF_btm_runfast_isoevents = cell(1, size(spk_inx_neurons,2));
    dfOvF_btm_runslow_isoevents = cell(1, size(spk_inx_neurons,2));
    dfOvF_btm_stay_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);% average df/f of isolated event of each cell, averaged across events, neurons*frames
    dfOvF_btm_runfast_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);
    dfOvF_btm_runslow_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);
    %     dfOvF_spk_run_begin_neurons = zeros(size(spk_inx_neurons,2),30);
    %     dfOvF_spk_run_middle_neurons = zeros(size(spk_inx_neurons,2),30);
    %     dfOvF_spk_run_late_neurons = zeros(size(spk_inx_neurons,2),30);
    dfOvF_btm_stay_rand_isoevents = [];
    dfOvF_btm_runfast_rand_isoevents = [];
    dfOvF_btm_runslow_rand_isoevents = [];
    nevent = zeros(1,size(spk_inx_neurons,2));
    neventrun = zeros(1,size(spk_inx_neurons,2));
    %     isospk_run_begin = 0;
    %     isospk_run_middle = 0;
    %     isospk_run_late = 0;
    %     dfOvF_spk_run_begin = cell(1,size(spk_inx_neurons,2));
    %     dfOvF_spk_run_middle = cell(1,size(spk_inx_neurons,2));
    %     dfOvF_spk_run_late = cell(1,size(spk_inx_neurons,2));
    for c = 1: size(spk_inx_neurons,2)% for each cell
        %stationary
        spk_iso_stay_neurons{c} = intersect(stay_vec,isolated_inx_neurons{c});
        %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
        isoevent_stay_neurons{c} = zeros(length(spk_iso_stay_neurons{c}),30); % for each cell, create a matrix, # of events*30(1s), events*frames, each row is a stationary event
        dfOvF_btm_stay_isoevents{c} =  zeros(length(spk_iso_stay_neurons{c}),30);% number of events*frames
        dfOvF_c = dfOvF_btm(:,c);
        for s = 1:length(spk_iso_stay_neurons{c})
            isoevent_stay_neurons{c}(s,:) =  spk_iso_stay_neurons{c}(s)-11:spk_iso_stay_neurons{c}(s)+18;
            dfOvF_btm_stay_isoevents{c}(s,:) = dfOvF_c(isoevent_stay_neurons{c}(s,:));
        end
        
        %running
        spk_iso_runfast_neurons{c} = intersect(run_fast_vec,isolated_inx_neurons{c});
        %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
        isoevent_runfast_neurons{c} = zeros(length(spk_iso_runfast_neurons{c}),30);
        dfOvF_btm_runfast_isoevents{c} =  zeros(length(spk_iso_runfast_neurons{c}),30); %each row is a spike, each column is a frame
        for s = 1:length(spk_iso_runfast_neurons{c}) % for each spike of each cell
            isoevent_runfast_neurons{c}(s,:) =  spk_iso_runfast_neurons{c}(s)-11:spk_iso_runfast_neurons{c}(s)+18;
            dfOvF_btm_runfast_isoevents{c}(s,:) = dfOvF_c(isoevent_runfast_neurons{c}(s,:));
            %-------------------------------------------------------------------------------
            % decide when the spike happens during running and
            % catagorize them based on when it happens --- for comparing the amplitude of isolated events during different running period
            % first find the first frame of that running window,
%             a = frm_run_starts(frm_run_starts <= spk_iso_runfast_neurons{c}(s)); % these are incides of all first frames of running that is before the spike
%             tfromstart = spk_iso_runfast_neurons{c}(s) - a(end); %the last one in a is the first frame that is in the same running window as the spike
%             if tfromstart <= 15
%                 isospk_run_begin = isospk_run_begin + 1;
%                 dfOvF_spk_run_begin{c} = cat(1,dfOvF_spk_run_begin{c},dfOvF_btm_runfast_isoevents{c}(s,:));
%             elseif tfromstart > 15 && tfromstart <= 30
%                 isospk_run_middle = isospk_run_middle + 1;
%                 dfOvF_spk_run_middle{c} = cat(1,dfOvF_spk_run_middle{c},dfOvF_btm_runfast_isoevents{c}(s,:));
%             else
%                 isospk_run_late = isospk_run_late + 1;
%                 dfOvF_spk_run_late{c} = cat(1,dfOvF_spk_run_late{c},dfOvF_btm_runfast_isoevents{c}(s,:));
%             end
        end
        spk_iso_runslow_neurons{c} = intersect(run_slow_vec,isolated_inx_neurons{c});
        %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
        isoevent_runslow_neurons{c} = zeros(length(spk_iso_runslow_neurons{c}),30);
        dfOvF_btm_runslow_isoevents{c} =  zeros(length(spk_iso_runslow_neurons{c}),30); %each row is a spike, each column is a frame
        for s = 1:length(spk_iso_runslow_neurons{c}) % for each spike of each cell
            isoevent_runslow_neurons{c}(s,:) =  spk_iso_runslow_neurons{c}(s)-11:spk_iso_runslow_neurons{c}(s)+18;
            dfOvF_btm_runslow_isoevents{c}(s,:) = dfOvF_c(isoevent_runslow_neurons{c}(s,:));
        end
        
        % draw same number of events from stationary and running and then do average-----------------------------------------------------
        nevent_stay = size(dfOvF_btm_stay_isoevents{c},1);
        nevent_runfast = size(dfOvF_btm_runfast_isoevents{c},1);
        nevent_runslow = size(dfOvF_btm_runslow_isoevents{c},1);
        nevent(c) = min([nevent_stay,nevent_runfast,nevent_runslow]);
       % if nevent(c)>0
            rand_stay = randperm(nevent_stay,nevent(c));
            rand_runfast = randperm(nevent_runfast,nevent(c));
            rand_runslow = randperm(nevent_runslow,nevent(c));
            dfOvF_btm_stay_rand_isoevents = dfOvF_btm_stay_isoevents{c}(rand_stay,:);
            dfOvF_btm_stay_isoevents_neurons(c,:) = mean(dfOvF_btm_stay_rand_isoevents); % average across trials
            dfOvF_btm_runfast_rand_isoevents = dfOvF_btm_runfast_isoevents{c}(rand_runfast,:);
            dfOvF_btm_runfast_isoevents_neurons(c,:) = mean(dfOvF_btm_runfast_rand_isoevents);
            dfOvF_btm_runslow_rand_isoevents = dfOvF_btm_runslow_isoevents{c}(rand_runslow,:);
            dfOvF_btm_runslow_isoevents_neurons(c,:) = mean(dfOvF_btm_runslow_rand_isoevents);
       % end
        
%         %for running, draw same number of events from different running times and then do average
%         if isempty(dfOvF_spk_run_begin{c}) == 0 && isempty(dfOvF_spk_run_middle{c}) == 0 && isempty(dfOvF_spk_run_late{c}) == 0 % for each cell, if there is isolated spike during every running period
%             nevent_begin = size(dfOvF_spk_run_begin{c},1);
%             nevent_middle = size(dfOvF_spk_run_middle{c},1);
%             nevent_late = size(dfOvF_spk_run_late{c},1);
%             neventrun(c) = min([nevent_begin,nevent_middle,nevent_late]);
%             rand_begin = randperm(nevent_begin,neventrun(c));
%             rand_middle = randperm(nevent_middle,neventrun(c));
%             rand_late = randperm(nevent_late,neventrun(c));
%             dfOvF_spk_run_begin_rand = dfOvF_spk_run_begin{c}(rand_begin,:);
%             dfOvF_spk_run_begin_neurons(c,:) = mean(dfOvF_spk_run_begin_rand); % average across trials
%             dfOvF_spk_run_middle_rand = dfOvF_spk_run_middle{c}(rand_middle,:);
%             dfOvF_spk_run_middle_neurons(c,:) = mean(dfOvF_spk_run_middle_rand);
%             dfOvF_spk_run_late_rand = dfOvF_spk_run_late{c}(rand_late,:);
%             dfOvF_spk_run_late_neurons(c,:) = mean(dfOvF_spk_run_late_rand);
%         end
        
    end
    
    % if there's a cell that doesn't spike during running,there will be a line of NaNs in dfOvF_btm_run/stay_isoevents_neurons, need to delete this Nan line before doing average
    dfOvF_btm_runfast_isoevents_neurons = dfOvF_btm_runfast_isoevents_neurons(all(~isnan(dfOvF_btm_runfast_isoevents_neurons),2),:);
    ave_dfOvF_btm_runfast_iso_session = mean(dfOvF_btm_runfast_isoevents_neurons);% average across cells
    ste_dfOvF_btm_runfast_iso_session = std(dfOvF_btm_runfast_isoevents_neurons,0,1)/sqrt(size(dfOvF_btm_runfast_isoevents_neurons,1)); % ste
    
    dfOvF_btm_runslow_isoevents_neurons = dfOvF_btm_runslow_isoevents_neurons(all(~isnan(dfOvF_btm_runslow_isoevents_neurons),2),:);
    ave_dfOvF_btm_runslow_iso_session = mean(dfOvF_btm_runslow_isoevents_neurons);% average across cells
    ste_dfOvF_btm_runslow_iso_session = std(dfOvF_btm_runslow_isoevents_neurons,0,1)/sqrt(size(dfOvF_btm_runslow_isoevents_neurons,1));
    
    dfOvF_btm_stay_isoevents_neurons = dfOvF_btm_stay_isoevents_neurons(all(~isnan(dfOvF_btm_stay_isoevents_neurons),2),:);
    ave_dfOvF_btm_stay_iso_session = mean(dfOvF_btm_stay_isoevents_neurons);
    ste_dfOvF_btm_stay_iso_session = std(dfOvF_btm_stay_isoevents_neurons,0,1)/sqrt(size(dfOvF_btm_stay_isoevents_neurons,1));
    
%     % for dfOvF_spk_run_begin/middle/late_neurons, delete those lines with zeros, that means that cell doesn't fire during all 3 conditions
%     dfOvF_spk_run_begin_neurons = dfOvF_spk_run_begin_neurons(any(dfOvF_spk_run_begin_neurons~=0,2),:);
%     dfOvF_spk_run_middle_neurons = dfOvF_spk_run_middle_neurons(any(dfOvF_spk_run_middle_neurons~=0,2),:);
%     dfOvF_spk_run_late_neurons = dfOvF_spk_run_late_neurons(any(dfOvF_spk_run_late_neurons~=0,2),:);
%     
    %     %count how many events you have in total
    %     nevents_run = sum(cellfun(@(x)x(1),cellfun(@size,isoevent_run_neurons,'UniformOutput',false)));% # of isolated events for all cells during running
    %     nevents_stay = sum(cellfun(@(x)x(1),cellfun(@size,isoevent_stay_neurons,'UniformOutput',false)));
    %
    %plot
    x = (1:30)/30;
    Caevent = figure;
    errorbar(x,ave_dfOvF_btm_stay_iso_session,ste_dfOvF_btm_stay_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
    errorbar(x,ave_dfOvF_btm_runfast_iso_session,ste_dfOvF_btm_runfast_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
    errorbar(x,ave_dfOvF_btm_runslow_iso_session,ste_dfOvF_btm_runslow_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;

    % number of calcium events is wrong in the current legend
    legend('stationary', 'fast running','slow running');
    ylabel('df/f');
    xlabel('time(s)');
    title(sessions(ii));
    text(0.75,0.45, ['nevents = ' num2str(sum(nevent))]);
    text(0.75,0.4, ['nPCs = ' num2str(size(dfOvF_btm_stay_isoevents_neurons,1))]);
    saveas(Caevent,[image_analysis_dest '\' days{ii} '_CaAmp']);
    % save variables
    
%     Caevent_run = figure;
%     fast_errbar(x,dfOvF_spk_run_begin_neurons,1,'color',[0.1373 0.5451 0.2706]);hold on;
%     fast_errbar(x,dfOvF_spk_run_middle_neurons,1,'color',[0.2549 0.6706 0.3647]); hold on;
%     fast_errbar(x,dfOvF_spk_run_late_neurons,1,'color',[0.7294 0.8941 0.7020]); hold on;
%     legend('0-0.5s','0.5-1s','later than 1s');
%     ylabel('df/f');
%     xlabel('time(s)');
%     title(sessions(ii));
%     text(0.75,0.4, ['nPCs = ' num2str(size(dfOvF_spk_run_begin_neurons,1))]);
%     saveas(Caevent_run,[image_analysis_dest '\' days{ii} '_CaAmp_run']);
%     
    % if isempty(frm_run_midpart_cell) == 0
    save([image_analysis_dest sessions{ii} '_isoCaEvent.mat'],'isolated_inx_neurons','spk_iso_stay_neurons','spk_iso_runfast_neurons',...
        'spk_iso_runslow_neurons','isoevent_stay_neurons','isoevent_runfast_neurons','isoevent_runslow_neurons','dfOvF_btm_stay_isoevents',...
        'dfOvF_btm_runfast_isoevents','dfOvF_btm_runslow_isoevents',...
        'dfOvF_btm_stay_isoevents_neurons','dfOvF_btm_runfast_isoevents_neurons',...
        'dfOvF_btm_runslow_isoevents_neurons','ave_dfOvF_btm_stay_iso_session','ste_dfOvF_btm_stay_iso_session',...
        'ave_dfOvF_btm_runfast_iso_session','ste_dfOvF_btm_runfast_iso_session',...
        'ave_dfOvF_btm_runslow_iso_session','ste_dfOvF_btm_runslow_iso_session','nevent');
    % else
    %     fprintf(['no running trials longer than 2s in this session ' sessions{ii}]);
    %     save([image_analysis_dest sessions{ii} '_isoCaEvent.mat'],'frm_stay_midpart',...
    %         'frm_run_midpart','isolated_inx_neurons','spk_iso_stay_neurons',...
    %         'isoevent_stay_neurons','dfOvF_btm_stay_isoevents','dfOvF_btm_stay_isoevents_neurons',...
    %         'ave_dfOvF_btm_stay_iso_session','ste_dfOvF_btm_stay_iso_session','nevents_stay');
end
    
% %%
% % decide when does the isolated spike happen during running? (beginning?middle?late?)
% % spk_iso_run_neurons tells you the frame index of the peak of isolated spike,
% % diff between spk_iso_run_neurons and the first frame of each running window tells you in which part it happens.
% % for each spk event, first find which running window it belongs to.
% % and is the amplitude influenced by when it happens during running?
% for   ii = 1: length(sessions)
%     image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
%     behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{ii} '\'];
%     behav_output = load([behav_dest days{ii} '_behavAnalysis.mat']);
%     frm_run_cell = behav_output.frames_run_cell;
%     frm_run_starts = zeros(1,size(frm_run_cell,2));
%     isoCa_output = load([image_analysis_dest sessions{ii} '_isoCaEvent.mat']);
%     spk_iso_run_neurons_vec = isoCa_output.spk_iso_run_neurons_vec;
%     for r = 1:size(frm_run_cell,2)
%         frm_run_starts(r) = frm_run_cell{r}(1);
%     end
%     isospk_run_begin = 0;
%     isospk_run_middle = 0;
%     isospk_run_late = 0;
%     for s = 1:length(spk_iso_run_neurons_vec) % for each spike index (this vector contains all indices of all spikes of all neurons)
%         % find the first frame of that running window,
%         a = frm_run_starts(frm_run_starts <= spk_iso_run_neurons_vec(s)); % these are incides of all first frames of running that is before the spike
%         tfromstart = spk_iso_run_neurons_vec(s) - a(end); %the last one in a is the first frame that is in the same running window as the spike
%         if tfromstart <= 15
%             isospk_run_begin = isospk_run_begin + 1;
%         elseif tfromstart > 15 && tfromstart <= 30
%             isospk_run_middle = isospk_run_middle + 1;
%         else
%             isospk_run_late = isospk_run_late + 1;
%         end
%     end
%     save([image_analysis_dest sessions{ii} '_isoCaEvent.mat'],'isospk_run_begin','isospk_run_middle','isospk_run_late','-append');
%     
% end
% 
%%
% decide when does the isolated spike happen during running? (beginning?middle?late?)
% spk_iso_run_neurons tells you the frame index of the peak of isolated spike,
% diff between spk_iso_run_neurons and the first frame of each running window tells you in which part it happens.
% for each spk event, first find which running window it belongs to. 

