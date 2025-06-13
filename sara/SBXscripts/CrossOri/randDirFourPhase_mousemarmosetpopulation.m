clc;clear all;close all
%% pull mouse data
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC','SCN'); 
area = 'all_areas';
area_list = strvcat('V1','V1');
narea = length(area_list);
nCells = [];


for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];

    pattern_ind = (Zp_all-Zc_all);    
    pattind = mean(pattern_ind,1);
    pattpeak = b_all+amp_all;
    ZpZcpeak = max(pattern_ind,[],1);

    nan_ind = isnan(amp_all);
    amp_all(nan_ind)=0;
    med_amp = median(amp_all(resp_ind));
    med_peak = median(pattpeak(resp_ind));

    figure(1)
        subplot(3,3,iA)
            scatter((-1*amp_all(resp_ind)),pattpeak(resp_ind),7,'filled')
            hold on
            yline(med_peak)
            xline(-1*med_amp)
            if iA <3
                ylim([-5 10])
                xlim([-6 0])
            else; end
            xlabel('Spatial invariance (-Amplitude)')
            ylabel('Direction invariance (Fit peak)')
            set(gca,'TickDir','out'); axis square
            if iA ==1;  subtitle(['mouse V1 L2/3']); else; subtitle(['mouse V1 L4']); end
         subplot(3,3,6)
            scatter((-1*amp_all(resp_ind)),pattpeak(resp_ind),7,'filled')
            hold on
            yline(med_peak)
            xline(-1*med_amp)
            if iA <3
                ylim([-5 10])
                xlim([-6 0])
            else; end
            if iA >1
                yline(med_peak,'r')
                xline(-1*med_amp,'r')
            else; end
            xlabel('Spatial invariance (-Amplitude)')
            ylabel('Direction invariance (Fit peak)')
            set(gca,'TickDir','out'); axis square
            subtitle(['mouse V1 - L2/3-blue, L4-red'])
       if iA ==1
         subplot(3,3,4)
            amp_twothree = amp_all;
            peak_twothree=pattpeak;
            resp_twothree = resp_ind;
            scatter((-1*amp_all(resp_ind)),pattpeak(resp_ind),7,'filled')
            hold on
                yline(med_peak)
                xline(-1*med_amp)    
                ylim([-5 10])
                xlim([-6 0])
            subtitle(['mouse L2/3-blue, marm-red'])
            xlabel('Spatial invariance (-Amplitude)');ylabel('Direction invariance (Fit peak)');set(gca,'TickDir','out'); axis square
       end
        if iA ==2
         subplot(3,3,5)
            scatter((-1*amp_all(resp_ind)),pattpeak(resp_ind),7,'filled')
            hold on
                yline(med_peak)
                xline(-1*med_amp)
                ylim([-5 10])
                xlim([-6 0])
            subtitle(['mouse L4-blue, marm-red'])
            xlabel('Spatial invariance (-Amplitude)');ylabel('Direction invariance (Fit peak)');set(gca,'TickDir','out'); axis square
        end
if iA == 1
    figure(2)
        subplot(4,4,1)
        scatter(Zc_all(1,resp_ind),Zp_all(1,resp_ind),7,'filled')
        hold on
else; end
end
    % amp_layers = [amp_all;amp_twothree];
    % peak_layers = [pattpeak;peak_twothree];
    % resp_layers = [resp_ind; length(amp_all)+resp_twothree];
    % med_al = median(amp_layers(resp_layers));
    % med_pl = median(peak_layers(resp_layers));
    % figure(1);
    % subplot(3,3,6)
    %     scatter((-1*(amp_layers(resp_layers))),peak_layers(resp_layers),7,'filled')
    %     hold on
    %     yline(med_pl)
    %     xline(-1*med_al)
    %     xlabel('Spatial invariance (-Amplitude)')
    %         ylabel('Direction invariance (Fit peak)')
    %         set(gca,'TickDir','out'); axis square

%% load marm data
clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = ([base '\Analysis\Neuropixel\CrossOri\randDirFourPhase\summaries']);
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randDirFourPhase_CrossOri';

load(fullfile(summaryDir,[svName '_Summary_V1_MAR.mat']))

    pattern_ind = (Zp_all-Zc_all);    
    pattind = mean(pattern_ind,1);
    pattpeak = max(pattern_ind,[],1);   

    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));

    nan_ind = isnan(amp_all);
    amp_all(nan_ind)=0;
    med_amp = median(amp_all(resp_ind));
    med_peak = median(pattpeak(resp_ind));

figure(1)
    subplot(3,3,3)
        scatter((-1*(amp_all(resp_ind))),pattpeak(resp_ind),7,'filled')
        hold on
        xlim([-6 0])
        ylim([-5 10])
        xlabel('Spatial invariance (-Amplitude)')
        ylabel('Direction invariance (Fit peak)')
        set(gca,'TickDir','out'); axis square
        subtitle(['marmoset V1'])
    subplot(3,3,4)
        scatter((-1*(amp_all(resp_ind))),pattpeak(resp_ind),7,'filled')
        hold on
        yline(med_peak,'r')
        xline(-1*med_amp,'r')
        xlim([-6 0])
        ylim([-5 10])
    subplot(3,3,5)
        scatter((-1*(amp_all(resp_ind))),pattpeak(resp_ind),7,'filled')
        hold on
        yline(med_peak,'r')
        xline(-1*med_amp,'r')
        xlim([-6 0])
        ylim([-5 10])
    % subplot(3,3,6)
    %     scatter((-1*(amp_all(resp_ind))),pattpeak(resp_ind),7,'filled')
    %     hold on
    %     yline(med_peak,'r')
    %     xline(-1*med_amp,'r')
    %     xlim([-6 0])
    %     ylim([-5 10])
figure(2)
    subplot(4,4,1)
        scatter(Zc_all(4,resp_ind),Zp_all(4,resp_ind),7,'filled')
        hold on
        xlabel('Zc'); ylabel('Zp'); xlim([-4 8]); ylim([-4 8]); plotZcZpBorders; set(gca,'TickDir','out'); axis square


%% print figs
stop
    figure(1); print(fullfile(outDir, [svName '_mousemarmoset_population.pdf']),'-dpdf', '-fillpage') 
    figure(2); print(fullfile(outDir, [svName '_mousemarmoset_populationZpZc.pdf']),'-dpdf', '-fillpage') 


%%

figure;
scatter(DSI_all,b_all)
