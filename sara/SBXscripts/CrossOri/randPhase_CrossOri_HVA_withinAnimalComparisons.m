close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Figures', 'Lab_meeting','250218');
svName = 'randPhase';
driver = strvcat('SLC','SLC','SLC','SLC'); 
area = 'all_areas';
area_list = strvcat('V1','LM','AL','PM');
narea = length(area_list);
nCells = [];



%% i1380

figure(1);
for iA = 1
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    fprintf([area_list(iA,:) '\n'])
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end

   ie = 3; % 230728_i1380
        sig_stim    = cell2mat(byexpt(ie).sig_stim);
        sig_dir     = cell2mat(byexpt(ie).sig_dir);
        DSI         = find(cell2mat(byexpt(ie).DSI)>0.5);
        ind_expt    = intersect(intersect(sig_stim,sig_dir),DSI);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,1)
            cdfplot(amp(ind_expt));
            hold on
            % xlabel('fit amp1');


    load(fullfile(summaryDir,[svName '_OnePhaseSummary_' area_list(iA,:) '_' driver(iA,:) '.mat']))
    fprintf([area_list(iA,:) ' ncells=' num2str(length(ind_all)) '\n'])
    
    ie = 5; % 240109_i1380
        ind_expt    = cell2mat(byexpt(ie).ind);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,1)
            cdfplot(amp(ind_expt));
            hold on
            % xlabel('fit amp1'); 
            subtitle([area_list(iA,:) ' i1380'])
            xlim([0 5]); set(gca,'TickDir','out'); box off; axis square
end


for iA = 2
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    fprintf([area_list(iA,:) '\n'])
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end

   ie = 2; % 230913_i1380
        sig_stim    = cell2mat(byexpt(ie).sig_stim);
        sig_dir     = cell2mat(byexpt(ie).sig_dir);
        DSI         = find(cell2mat(byexpt(ie).DSI)>0.5);
        ind_expt    = intersect(intersect(sig_stim,sig_dir),DSI);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,2)
            cdfplot(amp(ind_expt));
            hold on
            % xlabel('fit amp1');


    load(fullfile(summaryDir,[svName '_OnePhaseSummary_' area_list(iA,:) '_' driver(iA,:) '.mat']))
    fprintf([area_list(iA,:) ' ncells=' num2str(length(ind_all)) '\n'])
    
    ie = 2; % 240116_i1380
        ind_expt    = cell2mat(byexpt(ie).ind);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,2)
            cdfplot(amp(ind_expt));
            hold on
            xlim([0 5]); set(gca,'TickDir','out'); box off; axis square
            % xlabel('fit amp1'); 
            subtitle([area_list(iA,:) ' i1380'])
end



%% i1396

figure(1);
for iA = 1
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    fprintf([area_list(iA,:) '\n'])
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end

   ie = 9; % 240507_i1396
        sig_stim    = cell2mat(byexpt(ie).sig_stim);
        sig_dir     = cell2mat(byexpt(ie).sig_dir);
        DSI         = find(cell2mat(byexpt(ie).DSI)>0.5);
        ind_expt    = intersect(intersect(sig_stim,sig_dir),DSI);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,3)
            cdfplot(amp(ind_expt));
            hold on
            xlim([0 5]); set(gca,'TickDir','out'); box off; axis square
            % xlabel('fit amp1');
            subtitle([area_list(iA,:) ' i1396'])

end


for iA = 2
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    fprintf([area_list(iA,:) '\n'])
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end

   ie = 4; % 240501_i1396
        sig_stim    = cell2mat(byexpt(ie).sig_stim);
        sig_dir     = cell2mat(byexpt(ie).sig_dir);
        DSI         = find(cell2mat(byexpt(ie).DSI)>0.5);
        ind_expt    = intersect(intersect(sig_stim,sig_dir),DSI);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,4)
            cdfplot(amp(ind_expt));
            hold on
            % xlabel('fit amp1');


    load(fullfile(summaryDir,[svName '_OnePhaseSummary_' area_list(iA,:) '_' driver(iA,:) '.mat']))
    fprintf([area_list(iA,:) ' ncells=' num2str(length(ind_all)) '\n'])
    
    ie = 4; % 240510_i1396
        ind_expt    = cell2mat(byexpt(ie).ind);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,4)
            cdfplot(amp(ind_expt));
            hold on
            xlim([0 5]); set(gca,'TickDir','out'); box off; axis square
            % xlabel('fit amp1'); 
            subtitle([area_list(iA,:) ' i1396'])
end

for iA = 3
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    fprintf([area_list(iA,:) '\n'])
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end

   ie = 5; % 240502_i1396
        sig_stim    = cell2mat(byexpt(ie).sig_stim);
        sig_dir     = cell2mat(byexpt(ie).sig_dir);
        DSI         = find(cell2mat(byexpt(ie).DSI)>0.5);
        ind_expt    = intersect(intersect(sig_stim,sig_dir),DSI);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,5)
            cdfplot(amp(ind_expt));
            hold on
            % xlabel('fit amp1');


    load(fullfile(summaryDir,[svName '_OnePhaseSummary_' area_list(iA,:) '_' driver(iA,:) '.mat']))
    fprintf([area_list(iA,:) ' ncells=' num2str(length(ind_all)) '\n'])
    
    ie = 1; % 240508_i1396
        ind_expt    = cell2mat(byexpt(ie).ind);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,5)
            cdfplot(amp(ind_expt));
            hold on
            xlim([0 5]); set(gca,'TickDir','out'); box off; axis square
            % xlabel('fit amp1'); 
            subtitle([area_list(iA,:) 'i1396'])
end



%%

figure(1);
for iA = 3
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    fprintf([area_list(iA,:) '\n'])
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end

   ie = 4; % 240409_i1394
        sig_stim    = cell2mat(byexpt(ie).sig_stim);
        sig_dir     = cell2mat(byexpt(ie).sig_dir);
        DSI         = find(cell2mat(byexpt(ie).DSI)>0.5);
        ind_expt    = intersect(intersect(sig_stim,sig_dir),DSI);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,6)
            cdfplot(amp(ind_expt));
            hold on
            % xlabel('fit amp1');

    ie = 6; % 240611_i1394
        sig_stim    = cell2mat(byexpt(ie).sig_stim);
        sig_dir     = cell2mat(byexpt(ie).sig_dir);
        DSI         = find(cell2mat(byexpt(ie).DSI)>0.5);
        ind_expt    = intersect(intersect(sig_stim,sig_dir),DSI);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,6)
            cdfplot(amp(ind_expt));
            hold on
            % xlabel('fit amp1');

    load(fullfile(summaryDir,[svName '_OnePhaseSummary_' area_list(iA,:) '_' driver(iA,:) '.mat']))
    fprintf([area_list(iA,:) ' ncells=' num2str(length(ind_all)) '\n'])
    
    ie = 3; % 240520_i1394
        ind_expt    = cell2mat(byexpt(ie).ind);
        amp         = cell2mat(byexpt(ie).amp);

        subplot(2,4,6)
            cdfplot(amp(ind_expt));
            hold on
            xlim([0 5]); set(gca,'TickDir','out'); box off; axis square
            % xlabel('fit amp1'); 
            subtitle([area_list(iA,:) 'i1394'])
end












%%
figure(1);  print(fullfile(outDir, [svName '_OnetoFourPhaseComparison_WithinArea_ByAnimal.pdf']),'-dpdf', '-fillpage') 