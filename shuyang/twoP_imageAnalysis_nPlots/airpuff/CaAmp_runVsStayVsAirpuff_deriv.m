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
sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042'};
days = {'1040-191114_1','1039-191115_1','1041-191115_1','1042-191115_1'};
%sessionID = {'1023-190923_1','','','','','','','',''};
image_analysis_base = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\'; 
color_code = {'c','r','y','g'};

%% SectionII: for each session: running experiment, event amplitude during running vs stationary
% get events from each session, and save all of them into one variable at
% the end of the for loop
for ii = 2: length(sessions)
   % load data
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\']; 
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvF_btm = dfOvF_strct.dfOvF_btm_cl;
    behav_dest = ['Z:\Analysis\Airpuff_analysis\behavioral_analysis\' days{ii} '\'];
    behav_output = load([behav_dest days{ii} '_behavAnalysis.mat']);
    %behav_output = load([behav_dest days{i} '_first18000frames_behavAnalysis.mat']);
    frm_stay = behav_output.frm_stay;
    airpuff_period = behav_output.airpuff_period;
    run_noairpuff = behav_output.run_noairpuff;
    stay_noairpuff = behav_output.stay_noairpuff;
    puffresp_stay = intersect(airpuff_period,frm_stay);
    
    % for each cell,find isolated events in running and stationary (no event within 600ms before that event)
    spk_deconv_output = load([image_analysis_dest '\derivative\' sessions{ii},'_spikes.mat']);
    spk_inx_neurons = spk_deconv_output.spk_inx;
    
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
    spk_iso_airpuff_stay_neurons = cell(1,size(spk_inx_neurons,2));
    isoevent_stay_neurons = cell(1, size(spk_inx_neurons,2));% index of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
    isoevent_run_neurons = cell(1, size(spk_inx_neurons,2));
    isoevent_airpuff_stay_neurons = cell(1, size(spk_inx_neurons,2));
    dfOvF_btm_stay_isoevents = cell(1, size(spk_inx_neurons,2)); %df/f of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
    dfOvF_btm_run_isoevents = cell(1, size(spk_inx_neurons,2));
    dfOvF_btm_airpuff_stay_isoevents = cell(1, size(spk_inx_neurons,2));
    dfOvF_btm_stay_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);% average df/f of isolated event of each cell, averaged across events, neurons*frames
    dfOvF_btm_run_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);
    dfOvF_btm_airpuff_stay_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);
    
    dfOvF_btm_stay_rand_isoevents = [];
    dfOvF_btm_run_rand_isoevents = [];
    nevent_stay = [];
    nevent_run = [];
    nevent = zeros(1,size(spk_inx_neurons,2));
    % write spike event into matrices after getting index of spike peak, each event is 30 frames
    for c = 1: size(spk_inx_neurons,2)
        %stationary
        spk_iso_stay_neurons{c} = intersect(stay_noairpuff,isolated_inx_neurons{c});
        %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
        isoevent_stay_neurons{c} = zeros(length(spk_iso_stay_neurons{c}),30); % for each cell, create a matrix, # of events*30(1s), events*frames, each row is a stationary event
        dfOvF_btm_stay_isoevents{c} =  zeros(length(spk_iso_stay_neurons{c}),30);% number of events*frames
        dfOvF_c = dfOvF_btm(:,c);
        for s = 1:length(spk_iso_stay_neurons{c})
             isoevent_stay_neurons{c}(s,:) =  spk_iso_stay_neurons{c}(s)-11:spk_iso_stay_neurons{c}(s)+18;
             dfOvF_btm_stay_isoevents{c}(s,:) = dfOvF_c(isoevent_stay_neurons{c}(s,:));
        end
        
        %running
        %if isempty(frm_run_midpart_cell) == 0 % if there are running trials that are longer than 2s, then continue
        spk_iso_run_neurons{c} = intersect(run_noairpuff,isolated_inx_neurons{c});
        %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
        isoevent_run_neurons{c} = zeros(length(spk_iso_run_neurons{c}),30);
        dfOvF_btm_run_isoevents{c} =  zeros(length(spk_iso_run_neurons{c}),30); %each row is a spike, each column is a frame
        for s = 1:length(spk_iso_run_neurons{c})
            isoevent_run_neurons{c}(s,:) =  spk_iso_run_neurons{c}(s)-11:spk_iso_run_neurons{c}(s)+18;
            dfOvF_btm_run_isoevents{c}(s,:) = dfOvF_c(isoevent_run_neurons{c}(s,:));
        end
        %airpuff
        spk_iso_airpuff_stay_neurons{c} = intersect(puffresp_stay,isolated_inx_neurons{c});
        %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
        isoevent_airpuff_stay_neurons{c} = zeros(length(spk_iso_airpuff_stay_neurons{c}),30);
        dfOvF_btm_airpuff_stay_isoevents{c} =  zeros(length(spk_iso_airpuff_stay_neurons{c}),30); %each row is a spike, each column is a frame
        for s = 1:length(spk_iso_airpuff_stay_neurons{c})
            isoevent_airpuff_stay_neurons{c}(s,:) =  spk_iso_airpuff_stay_neurons{c}(s)-11:spk_iso_airpuff_stay_neurons{c}(s)+18;
            dfOvF_btm_airpuff_stay_isoevents{c}(s,:) = dfOvF_c(isoevent_airpuff_stay_neurons{c}(s,:));
        end
        
        % draw same number of events from stationary and running and airpuff response then do average
        if isempty(dfOvF_btm_run_isoevents) == 0 % for each cell, if there is isolated spike during running
            nevent_stay = size(dfOvF_btm_stay_isoevents{c},1);
            nevent_run = size(dfOvF_btm_run_isoevents{c},1);
            nevent_puff_stay = size(dfOvF_btm_airpuff_stay_isoevents{c},1);
            nevent(c) = min([nevent_stay,nevent_run,nevent_puff_stay]);
            rand_stay = randperm(nevent_stay,nevent(c));
            rand_run = randperm(nevent_run,nevent(c));
            rand_puff_stay = randperm(nevent_puff_stay,nevent(c));
            dfOvF_btm_stay_rand_isoevents = dfOvF_btm_stay_isoevents{c}(rand_stay,:);     
            dfOvF_btm_stay_isoevents_neurons(c,:) = mean(dfOvF_btm_stay_rand_isoevents); % average across trials
            dfOvF_btm_run_rand_isoevents = dfOvF_btm_run_isoevents{c}(rand_run,:);
            dfOvF_btm_run_isoevents_neurons(c,:) = mean(dfOvF_btm_run_rand_isoevents);
            dfOvF_btm_puff_stay_rand_isoevents = dfOvF_btm_airpuff_stay_isoevents{c}(rand_puff_stay,:);
            dfOvF_btm_airpuff_stay_isoevents_neurons(c,:) = mean(dfOvF_btm_puff_stay_rand_isoevents);
        end
    end
           
    % if there's a cell that doesn't spike during running,there will be a line of NaNs in dfOvF_btm_run/stay_isoevents_neurons, need to delete this Nan line before doing average
    dfOvF_btm_run_isoevents_neurons = dfOvF_btm_run_isoevents_neurons(all(~isnan(dfOvF_btm_run_isoevents_neurons),2),:);
    ave_dfOvF_btm_run_iso_session = mean(dfOvF_btm_run_isoevents_neurons);% average across cells
    ste_dfOvF_btm_run_iso_session = std(dfOvF_btm_run_isoevents_neurons,0,1)/sqrt(size(dfOvF_btm_run_isoevents_neurons,1)); % ste
    
    dfOvF_btm_stay_isoevents_neurons = dfOvF_btm_stay_isoevents_neurons(all(~isnan(dfOvF_btm_stay_isoevents_neurons),2),:);
    ave_dfOvF_btm_stay_iso_session = mean(dfOvF_btm_stay_isoevents_neurons);
    ste_dfOvF_btm_stay_iso_session = std(dfOvF_btm_stay_isoevents_neurons,0,1)/sqrt(size(dfOvF_btm_stay_isoevents_neurons,1));
    
    
    dfOvF_btm_airpuff_stay_isoevents_neurons = dfOvF_btm_airpuff_stay_isoevents_neurons(all(~isnan(dfOvF_btm_airpuff_stay_isoevents_neurons),2),:);
    ave_dfOvF_btm_airpuff_stay_iso_session = mean(dfOvF_btm_airpuff_stay_isoevents_neurons);
    ste_dfOvF_btm_airpuff_stay_iso_session = std(dfOvF_btm_airpuff_stay_isoevents_neurons,0,1)/sqrt(size(dfOvF_btm_airpuff_stay_isoevents_neurons,1));
    
%     %count how many events you have in total
%     nevents_run = sum(cellfun(@(x)x(1),cellfun(@size,isoevent_run_neurons,'UniformOutput',false)));% # of isolated events for all cells during running
%     nevents_stay = sum(cellfun(@(x)x(1),cellfun(@size,isoevent_stay_neurons,'UniformOutput',false)));
%     
    %plot
    x = (1:30)/30;
    Caevent = figure;
    errorbar(x,ave_dfOvF_btm_stay_iso_session,ste_dfOvF_btm_stay_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
    errorbar(x,ave_dfOvF_btm_run_iso_session,ste_dfOvF_btm_run_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
    errorbar(x,ave_dfOvF_btm_airpuff_stay_iso_session,ste_dfOvF_btm_airpuff_stay_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20);hold on;
    % number of calcium events is wrong in the current legend
    legend('stationary', 'running','airpuff response_stationary');
    ylabel('df/f');
    xlabel('time(s)');
    title(sessions(ii));
    text(0.75,0.45, ['nevents = ' num2str(sum(nevent))]);
    text(0.75,0.4, ['nPCs = ' num2str(size(dfOvF_btm_stay_isoevents_neurons,1))]);
    saveas(Caevent,[image_analysis_dest '\derivative\' days{ii} '_CaAmp_runVsStayVsAirpuff']);
    % save variables
    
    % if isempty(frm_run_midpart_cell) == 0
    save([image_analysis_dest '\derivative\' sessions{ii} '_isoCaEvent.mat'],...
        'isolated_inx_neurons','spk_iso_stay_neurons','spk_iso_run_neurons',...
        'isoevent_stay_neurons','isoevent_run_neurons','dfOvF_btm_stay_isoevents','dfOvF_btm_run_isoevents',...
        'dfOvF_btm_stay_isoevents_neurons','dfOvF_btm_run_isoevents_neurons',...
        'ave_dfOvF_btm_stay_iso_session','ste_dfOvF_btm_stay_iso_session',...
        'ave_dfOvF_btm_run_iso_session','ste_dfOvF_btm_run_iso_session','nevent',...
        'dfOvF_btm_airpuff_stay_isoevents_neurons','ave_dfOvF_btm_airpuff_stay_iso_session','ste_dfOvF_btm_airpuff_stay_iso_session',...
        'dfOvF_btm_airpuff_stay_isoevents');
    % else
%     fprintf(['no running trials longer than 2s in this session ' sessions{ii}]);
%     save([image_analysis_dest sessions{ii} '_isoCaEvent.mat'],'frm_stay_midpart',...
%         'frm_run_midpart','isolated_inx_neurons','spk_iso_stay_neurons',...
%         'isoevent_stay_neurons','dfOvF_btm_stay_isoevents','dfOvF_btm_stay_isoevents_neurons',...
%         'ave_dfOvF_btm_stay_iso_session','ste_dfOvF_btm_stay_iso_session','nevents_stay');
    % end
end  

