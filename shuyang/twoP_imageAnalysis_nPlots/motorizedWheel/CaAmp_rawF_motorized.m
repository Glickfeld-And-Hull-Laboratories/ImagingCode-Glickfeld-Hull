% first part of this script is the same thing as CaAmp_runVsStay_motorized
% how is the Ca event amplitude different during the
% beginning of running, middle of running, and late running?
% also how many events are there during each running period

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
    rawF_output = load([image_analysis_dest sessions{ii},'_deconvolution_thresh-4_TCave_cl.mat']);
    rawF = rawF_output.TCave_cl;
    behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\running\behavioral_analysis\' days{ii} '\'];
    behav_output = load([behav_dest days{ii} '_behavAnalysis.mat']);
    run_fast_vec = behav_output.run_fast_vec;
    run_slow_vec = behav_output.run_slow_vec;
    stay_vec = behav_output.stay_vec;
    
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
    spk_iso_runfast_neurons = cell(1,size(spk_inx_neurons,2));
    spk_iso_runslow_neurons = cell(1,size(spk_inx_neurons,2));
    isoevent_stay_neurons = cell(1, size(spk_inx_neurons,2));% index of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
    isoevent_runfast_neurons = cell(1, size(spk_inx_neurons,2));
    isoevent_runslow_neurons = cell(1, size(spk_inx_neurons,2));
    rawF_stay_isoevents = cell(1, size(spk_inx_neurons,2)); %df/f of whole isolated event(1s) during stationary, 1*#of neurons, each cell element is a matrix(event*frames)
    rawF_runfast_isoevents = cell(1, size(spk_inx_neurons,2));
    rawF_runslow_isoevents = cell(1, size(spk_inx_neurons,2));
    rawF_stay_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);% average df/f of isolated event of each cell, averaged across events, neurons*frames
    rawF_runfast_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);
    rawF_runslow_isoevents_neurons = zeros(size(spk_inx_neurons,2),30);
  
    rawF_stay_rand_isoevents = [];
    rawF_runfast_rand_isoevents = [];
    rawF_runslow_rand_isoevents = [];
    nevent = zeros(1,size(spk_inx_neurons,2));
    neventrun = zeros(1,size(spk_inx_neurons,2));
  
    for c = 1: size(spk_inx_neurons,2)% for each cell
        %stationary
        spk_iso_stay_neurons{c} = intersect(stay_vec,isolated_inx_neurons{c});
        %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
        isoevent_stay_neurons{c} = zeros(length(spk_iso_stay_neurons{c}),30); % for each cell, create a matrix, # of events*30(1s), events*frames, each row is a stationary event
        rawF_stay_isoevents{c} =  zeros(length(spk_iso_stay_neurons{c}),30);% number of events*frames
        rawF_c = rawF(:,c);
        for s = 1:length(spk_iso_stay_neurons{c})
            isoevent_stay_neurons{c}(s,:) =  spk_iso_stay_neurons{c}(s)-11:spk_iso_stay_neurons{c}(s)+18;
            rawF_stay_isoevents{c}(s,:) = rawF_c(isoevent_stay_neurons{c}(s,:));
        end
        
        %running
        spk_iso_runfast_neurons{c} = intersect(run_fast_vec,isolated_inx_neurons{c});
        %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
        isoevent_runfast_neurons{c} = zeros(length(spk_iso_runfast_neurons{c}),30);
        rawF_runfast_isoevents{c} =  zeros(length(spk_iso_runfast_neurons{c}),30); %each row is a spike, each column is a frame
        for s = 1:length(spk_iso_runfast_neurons{c}) % for each spike of each cell
            isoevent_runfast_neurons{c}(s,:) =  spk_iso_runfast_neurons{c}(s)-11:spk_iso_runfast_neurons{c}(s)+18;
            rawF_runfast_isoevents{c}(s,:) = rawF_c(isoevent_runfast_neurons{c}(s,:));
        end
        spk_iso_runslow_neurons{c} = intersect(run_slow_vec,isolated_inx_neurons{c});
        %isolated event matrices, each event: 400ms (12 frames) before the peak and 600ms(18 frames) after the peak
        isoevent_runslow_neurons{c} = zeros(length(spk_iso_runslow_neurons{c}),30);
        rawF_runslow_isoevents{c} =  zeros(length(spk_iso_runslow_neurons{c}),30); %each row is a spike, each column is a frame
        for s = 1:length(spk_iso_runslow_neurons{c}) % for each spike of each cell
            isoevent_runslow_neurons{c}(s,:) =  spk_iso_runslow_neurons{c}(s)-11:spk_iso_runslow_neurons{c}(s)+18;
            rawF_runslow_isoevents{c}(s,:) = rawF_c(isoevent_runslow_neurons{c}(s,:));
        end
        
        % draw same number of events from stationary and running and then do average-----------------------------------------------------
        nevent_stay = size(rawF_stay_isoevents{c},1);
        nevent_runfast = size(rawF_runfast_isoevents{c},1);
        nevent_runslow = size(rawF_runslow_isoevents{c},1);
        nevent(c) = min([nevent_stay,nevent_runfast,nevent_runslow]);
       % if nevent(c)>0
            rand_stay = randperm(nevent_stay,nevent(c));
            rand_runfast = randperm(nevent_runfast,nevent(c));
            rand_runslow = randperm(nevent_runslow,nevent(c));
            rawF_stay_rand_isoevents = rawF_stay_isoevents{c}(rand_stay,:);
            rawF_stay_isoevents_neurons(c,:) = mean(rawF_stay_rand_isoevents); % average across trials
            rawF_runfast_rand_isoevents = rawF_runfast_isoevents{c}(rand_runfast,:);
            rawF_runfast_isoevents_neurons(c,:) = mean(rawF_runfast_rand_isoevents);
            rawF_runslow_rand_isoevents = rawF_runslow_isoevents{c}(rand_runslow,:);
            rawF_runslow_isoevents_neurons(c,:) = mean(rawF_runslow_rand_isoevents);
       % end
        
    end
    
    % if there's a cell that doesn't spike during running,there will be a line of NaNs in rawF_btm_run/stay_isoevents_neurons, need to delete this Nan line before doing average
    rawF_runfast_isoevents_neurons = rawF_runfast_isoevents_neurons(all(~isnan(rawF_runfast_isoevents_neurons),2),:);
    ave_rawF_runfast_iso_session = mean(rawF_runfast_isoevents_neurons);% average across cells
    ste_rawF_runfast_iso_session = std(rawF_runfast_isoevents_neurons,0,1)/sqrt(size(rawF_runfast_isoevents_neurons,1)); % ste
    
    rawF_runslow_isoevents_neurons = rawF_runslow_isoevents_neurons(all(~isnan(rawF_runslow_isoevents_neurons),2),:);
    ave_rawF_runslow_iso_session = mean(rawF_runslow_isoevents_neurons);% average across cells
    ste_rawF_runslow_iso_session = std(rawF_runslow_isoevents_neurons,0,1)/sqrt(size(rawF_runslow_isoevents_neurons,1));
    
    rawF_stay_isoevents_neurons = rawF_stay_isoevents_neurons(all(~isnan(rawF_stay_isoevents_neurons),2),:);
    ave_rawF_stay_iso_session = mean(rawF_stay_isoevents_neurons);
    ste_rawF_stay_iso_session = std(rawF_stay_isoevents_neurons,0,1)/sqrt(size(rawF_stay_isoevents_neurons,1));
    
    %plot
    x = (1:30)/30;
    Caevent = figure;
    errorbar(x,ave_rawF_stay_iso_session,ste_rawF_stay_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
    errorbar(x,ave_rawF_runfast_iso_session,ste_rawF_runfast_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;
    errorbar(x,ave_rawF_runslow_iso_session,ste_rawF_runslow_iso_session,'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20); hold on;

    % number of calcium events is wrong in the current legend
    legend('stationary', 'fast running','slow running');
    ylabel('df/f');
    xlabel('time(s)');
    title(sessions(ii));
    text(0.75,0.45, ['nevents = ' num2str(sum(nevent))]);
    text(0.75,0.4, ['nPCs = ' num2str(size(rawF_stay_isoevents_neurons,1))]);
    saveas(Caevent,[image_analysis_dest '\' days{ii} '_CaAmp_rawF']);

    save([image_analysis_dest sessions{ii} '_isoCaEvent.mat'],'rawF_stay_isoevents',...
        'rawF_runfast_isoevents','rawF_runslow_isoevents',...
        'rawF_stay_isoevents_neurons','rawF_runfast_isoevents_neurons',...
        'rawF_runslow_isoevents_neurons','ave_rawF_stay_iso_session','ste_rawF_stay_iso_session',...
        'ave_rawF_runfast_iso_session','ste_rawF_runfast_iso_session',...
        'ave_rawF_runslow_iso_session','ste_rawF_runslow_iso_session','nevent','-append');
   
end
  

%% calcium amplitude during different runing period

 for ii = 1: length(sessions)
     image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
     rawF_output = load([image_analysis_dest sessions{ii},'_deconvolution_thresh-4_TCave_cl.mat']);
     rawF = rawF_output.TCave_cl;
     behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\running\behavioral_analysis\' days{ii} '\'];
     behav_output = load([behav_dest days{ii} '_behavAnalysis.mat']);
     runstart = behav_output.run_start;
     Caoutput = load([image_analysis_dest sessions{ii} '_isoCaEvent.mat']);
     spk_iso_runslow_neurons = Caoutput.spk_iso_runslow_neurons;
     spk_iso_runfast_neurons = Caoutput.spk_iso_runfast_neurons;
     rawF_runfast_isoevents = Caoutput.rawF_runfast_isoevents;
     rawF_runslow_isoevents = Caoutput.rawF_runslow_isoevents;
     
     % put slow trials and fast trials together
     spk_iso_run_neurons = cell(size(spk_iso_runslow_neurons));
     rawF_run_isoevents = cell(size(spk_iso_runslow_neurons));
     for c = 1:size(spk_iso_runslow_neurons,2)
         spk_iso_run_neurons{c} = cat(1,spk_iso_runslow_neurons{c}, spk_iso_runfast_neurons{c}); 
         rawF_run_isoevents{c} = cat(1,rawF_runslow_isoevents{c},rawF_runfast_isoevents{c});
     end
     
     % decide when the spike happens during running and
     % catagorize them based on when it happens --- for comparing the amplitude of isolated events during different running period
     % first find the first frame of that running window

     rawF_spk_run_begin_neurons = zeros(size(spk_iso_run_neurons,2),30);
     rawF_spk_run_middle_neurons = zeros(size(spk_iso_run_neurons,2),30);
     rawF_spk_run_late_neurons = zeros(size(spk_iso_run_neurons,2),30);
     isospk_run_begin = 0;
     isospk_run_middle = 0;
     isospk_run_late = 0;
     rawF_spk_run_begin = cell(1,size(spk_iso_run_neurons,2));
     rawF_spk_run_middle = cell(1,size(spk_iso_run_neurons,2));
     rawF_spk_run_late = cell(1,size(spk_iso_run_neurons,2));
     for c = 1: size(spk_iso_runslow_neurons,2)
         for s = 1:length(spk_iso_run_neurons{c}) % for each spike of each cell
             a = runstart(runstart <= spk_iso_run_neurons{c}(s)); % these are incides of all first frames of running that is before the spike
             tfromstart = spk_iso_run_neurons{c}(s) - a(end); %the last one in a is the first frame that is in the same running window as the spike
             if tfromstart <= 15
                isospk_run_begin = isospk_run_begin + 1;
                rawF_spk_run_begin{c} = cat(1,rawF_spk_run_begin{c},rawF_run_isoevents{c}(s,:));
            elseif tfromstart > 15 && tfromstart <= 45
                isospk_run_middle = isospk_run_middle + 1;
                rawF_spk_run_middle{c} = cat(1,rawF_spk_run_middle{c},rawF_run_isoevents{c}(s,:));
            else
                isospk_run_late = isospk_run_late + 1;
                rawF_spk_run_late{c} = cat(1,rawF_spk_run_late{c},rawF_run_isoevents{c}(s,:));
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
     
     
     % for rawF_spk_run_begin/middle/late_neurons, delete those lines with zeros, that means that cell doesn't fire during all 3 conditions
     rawF_spk_run_begin_neurons = rawF_spk_run_begin_neurons(any(rawF_spk_run_begin_neurons~=0,2),:);
     rawF_spk_run_middle_neurons = rawF_spk_run_middle_neurons(any(rawF_spk_run_middle_neurons~=0,2),:);
     rawF_spk_run_late_neurons = rawF_spk_run_late_neurons(any(rawF_spk_run_late_neurons~=0,2),:);
     
     Caevent_run = figure;
     x = (1:30)/30;
     fast_errbar(x,rawF_spk_run_begin_neurons,1,'color',[0.1373 0.5451 0.2706]);hold on;
     fast_errbar(x,rawF_spk_run_middle_neurons,1,'color',[0.2549 0.6706 0.3647]); hold on;
     fast_errbar(x,rawF_spk_run_late_neurons,1,'color',[0.7294 0.8941 0.7020]); hold on;
     legend('0-0.5s','0.5-1.5s','later than 1.5s');
     ylabel('df/f');
     xlabel('time(s)');
     title(sessions(ii));
     text(0.75,0.4, ['nPCs = ' num2str(size(rawF_spk_run_begin_neurons,1))]);
     saveas(Caevent_run,[image_analysis_dest '\' days{ii} '_CaAmp_runPeriod_rawF']);
     
     Carun=figure;
     bar([isospk_run_begin,isospk_run_middle,isospk_run_late]);
     set(gca,'XTickLabel',{'0-0.5s','0.5-1.5s','later than 1.5s'});
     ylabel('number of events');
     saveas(Carun,[image_analysis_dest '\' days{ii} '_CaAmp_runPeriod_rawF']);

     save([image_analysis_dest sessions{ii} '_isoCaEvent.mat'],...
         'rawF_spk_run_begin_neurons','rawF_spk_run_middle_neurons','rawF_spk_run_late_neurons',...
         '-append');
     
 end