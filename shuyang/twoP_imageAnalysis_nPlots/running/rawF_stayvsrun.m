% compare the raw fluorescence during running vs. stationary WHEN THE CELL
% IS NOT FIRING. 
% doing this analysis because the baseline of CaAmp plot is different for
% stationary and running and I need to figure out why

clear;
sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; color_code = {'c','r','y','g'};

%% 
for ii = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    rawF_output = load([image_analysis_dest sessions{ii},'_deconvolution_thresh-4_TCave_cl.mat']);
    rawF = rawF_output.TCave_cl;
    CaAmp_output = load([image_analysis_dest sessions{ii} '_isoCaEvent.mat']);
    stay_CSpinx = CaAmp_output.spk_iso_stay_neurons; % when peak of Ca event happens
    run_CSpinx = CaAmp_output.spk_iso_run_neurons;
    
    stay_befoCS = cell(size(stay_CSpinx));
    run_befoCS = cell(size(stay_CSpinx));
    rawF_stay_befoCS = cell(size(stay_CSpinx));
    rawF_run_befoCS = cell(size(stay_CSpinx));
    rawF_stay_befoCS_vec = zeros(1,size(stay_CSpinx,2));
    rawF_run_befoCS_vec = zeros(1,size(stay_CSpinx,2));
    
    for c = 1:size(stay_CSpinx,2)
        for s1 = 1:length(stay_CSpinx{c}) %for each spike, take the baseline period before Ca transient
            stay_befoCS{c} = [stay_befoCS{c} stay_CSpinx{c}(s1)-18:stay_CSpinx{c}(s1)-4]; %peaks are at least 18 frames from each other (600ms, each frame is 33Hz.), and there're about 4 raising frames
        end
        
        for s2 = 1:length(run_CSpinx{c}) %for each spike, take the baseline period before Ca transient
            run_befoCS{c} = [run_befoCS{c} run_CSpinx{c}(s2)-18:run_CSpinx{c}(s2)-4];
        end
        
        rawF_stay_befoCS{c} = rawF(stay_befoCS{c},c);
        rawF_run_befoCS{c} = rawF(run_befoCS{c},c);
        
        rawF_stay_befoCS_vec(c) = mean(rawF_stay_befoCS{c}); %for each cell, average across time
        rawF_run_befoCS_vec(c) = mean(rawF_run_befoCS{c}); %for each cell, average across time
    end
    
    rawF_compare = [rawF_stay_befoCS_vec',rawF_run_befoCS_vec'];
    rawF_compare(~all(~isnan(rawF_compare),2),:)=[];
    %not all of the cells fire during slow running and fast running, so there will be lines in rawF_compare has NaN. delete those lines/cells
    
    
    x = [1,2];
    x_plot = repmat(x,size(rawF_compare,1),1);
    figure; 
    plot(x_plot',rawF_compare','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20,'color',[0.5294 0.5294 0.5294]); hold on;
    fast_errbar(x,rawF_compare,1,'color',[0.8896 0.1922 0.1557]);
    xlabel ('behavioral state');xlim([0.5 2.5]); axis square;
    set(gca,'XTick',x,'XTicklabel',{'stationary','running'});
    ylabel('raw F w/o CS');
    title(sessions{ii});
    savefig([image_analysis_dest sessions{ii} '_rawF_baseline.fig']);
    save([image_analysis_dest sessions{ii},'_deconvolution_thresh-4_TCave_cl.mat'],'rawF_compare','-append');
end


%% across sessions
rawF_compare_all = [];

for ii = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    rawF_output = load([image_analysis_dest sessions{ii},'_deconvolution_thresh-4_TCave_cl.mat']);
    rawF_compare = rawF_output.rawF_compare;
    rawF_compare_all = cat(1,rawF_compare_all,rawF_compare);
end

x = [1,2];
figure;
x_plot = repmat(x,size(rawF_compare_all,1),1);
plot(x_plot',rawF_compare_all','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20,'color',[0.5294 0.5294 0.5294]); hold on;
fast_errbar(x,rawF_compare_all,1,'color',[0.8896 0.1922 0.1557]);
xlabel ('behavioral state');xlim([0.5 2.5]); axis square;
set(gca,'XTick',x,'XTicklabel',{'stationary','running'});
ylabel('raw F without CS');
title('across sessions');
savefig(['Z:\Analysis\2P_MovingDots_Analysis\across_sessions\' 'across_sessions_rawF_baseline.fig']);



