
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
for ii = 4
   % load data
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    %load data, use data of only good cells
    rawF_output = load([image_analysis_dest sessions{ii}, '_deconvolution_thresh-4_TCave_cl.mat']);
    rawF = rawF_output.TCave_cl;
    behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days{ii} '\'];
    behav_output = load([behav_dest days{ii} '_behavAnalysis.mat']);
    %behav_output = load([behav_dest days{i} '_first18000frames_behavAnalysis.mat']);
    frm_stay_cell = behav_output.frames_stay_cell;
    frm_run_cell = behav_output.frames_run_cell;
    frm_stay = cell2mat(frm_stay_cell);
    frm_run = cell2mat(frm_run_cell);
    
    % for each cell,find isolated events in running and stationary (no event within 650ms before that event)
    spk_deriv_output = load([image_analysis_dest 'derivative\' sessions{ii}, '_spikes' '.mat']);
    spk_inx_cl = spk_deriv_output.spk_inx;
    peaks_run = [];
    peaks_stay = [];
    for c = 1: size(spk_inx_cl,2)
        spk_run = intersect(frm_run,spk_inx_cl{c});
        spk_stay = intersect(frm_stay,spk_inx_cl{c});
        peaks_run = cat(1,peaks_run,rawF(spk_run,c));
        peaks_stay = cat(1,peaks_stay,rawF(spk_stay,c));
    end
    
    hist_fig = figure;
    hist_fig.Units = 'centimeters';
    hist_fig.Position = [1 1 7.5 6];
    hold on;subplot(1,2,1);
    h1 = histogram(peaks_run,'BinWidth',200);
    h1.FaceColor = 'b'; %This way the color is the most different ,don't know how to change the colors respectively
    h1.EdgeColor = 'b';
    xlim([0 max(peaks_stay)]); 
    xlabel('rawF of peaks');
    ylabel('number of events');
    %title('df/F of spike peaks running');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',7);
    subplot(1,2,2);
    h2 = histogram(peaks_stay,'BinWidth',200);
    h2.FaceColor = 'r';
    h2.EdgeColor = 'r';
    xlim([0 max(peaks_stay)]); 
    title('stay');
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',7);
    hold off;
    saveas(hist_fig,[image_analysis_dest 'derivative\' sessions{ii} '_SpkeventAmp_hist']);
    fig_name = [sessions{ii} '_SpkeventAmp_hist_derivative'];
    path = [image_analysis_dest 'derivative\'];
    orient(hist_fig,'landscape');
    print(hist_fig,[path,fig_name],'-r600','-dpdf');

end
    


