close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randPhase';
driver = strvcat('Scn', 'SLC'); 
area = 'all_areas';
area_list = strvcat('V1', 'V1');
narea = length(area_list);
nCells = [];


figure(1);
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));
    if exist('red_cells_all','var')
        resp_ind = setdiff(resp_ind, red_cells_all);
    end
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(resp_ind))];



    nCellsResp = length(resp_ind);
    neuron_classes=[];
    clear ind_patt ind_comp

    for ip = 1:4
        ind_patt{ip} = intersect(find(Zp_all(ip,:)>1.28),find(Zp_all(ip,:)-Zc_all(ip,:)>1.28));
        ind_comp{ip} = intersect(find(Zc_all(ip,:)>1.28),find(Zc_all(ip,:)-Zp_all(ip,:)>1.28));
        ind_uncl{ip} = setdiff(1:length(Zc_all), [ind_patt{ip} ind_comp{ip}]);

        % Assign label for each neuron
        for ic = 1:nCellsResp
            i = resp_ind(ic);
            if ismember(i, ind_patt{ip})
                neuron_classes{ic}(ip) = 'p';
            elseif ismember(i, ind_comp{ip})
                neuron_classes{ic}(ip) = 'c';
            elseif ismember(i, ind_uncl{ip})
                neuron_classes{ic}(ip) = 'r';
            else
                neuron_classes{ic}(ip) = 'x'; % fallback in case itcs none
                error(['cell ' num2str(i) ' unclassed at phase=' num2str(ip)])
            end
        end
    end
   
    % Make lists for categories
    truly_p = [];
    truly_c = [];
    truly_r = [];
    switch_pc = [];
    switch_pr = [];
    switch_cr = [];
    switch_pcr = [];

    for i = 1:nCellsResp
        classes = unique(neuron_classes{i});
        classes(classes == 'x') = []; % remove unknowns if any
    
        switch length(classes)
            case 1
                switch classes
                    case 'p'
                        truly_p(end+1) = i;
                    case 'c'
                        truly_c(end+1) = i;
                    case 'r'
                        truly_r(end+1) = i;
                end
            case 2
                if all(ismember(classes, ['p' 'c']))
                    switch_pc(end+1) = i;
                elseif all(ismember(classes, ['p' 'r']))
                    switch_pr(end+1) = i;
                elseif all(ismember(classes, ['c' 'r']))
                    switch_cr(end+1) = i;
                end
            case 3
                switch_pcr(end+1) = i;
        end
    end

   % Labels and counts
    categories = {'Truly P', 'Truly C', 'Truly R', 'Switch P/C', 'Switch P/R', 'Switch C/R', 'Switch P/C/R'};
    counts = [length(truly_p), length(truly_c), length(truly_r), [length(switch_pc)+ length(switch_pr)+length(switch_cr)+length(switch_pcr)]];
    countsAll = [length(truly_p), length(truly_c), length(truly_r), ...
              length(switch_pc), length(switch_pr), length(switch_cr), length(switch_pcr)];
    
    % Convert to percentages
    percentages = 100 * counts / nCellsResp;
    
    % Plot
    figure(1);
    subplot(2,3,iA)
    bar(percentages);
    set(gca, 'XTickLabel', categories, 'XTickLabelRotation', 45);
    ylabel('percentage of cells');
    ylim([0 100]);
    subtitle(['mouse'  leg_str{iA}])
    movegui('center');
    hold on
    sgtitle('Cell classification switching: p=pattern, c=component, r=unclassified');




    %% pie charts


clear ind_patt ind_comp

    ind_patt(1) = length(intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28)),resp_ind));
    ind_patt(2) = length(intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),resp_ind));
    ind_patt(3) = length(intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(2,:)>1.28)),resp_ind));
    ind_patt(4) = length(intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),resp_ind));
    ind_comp(1) = length(intersect(intersect(find(Zc_all(1,:)>1.28),find(Zc_all(1,:)-Zp_all(1,:)>1.28)),resp_ind));
    ind_comp(2) = length(intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),resp_ind));
    ind_comp(3) = length(intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(3,:)-Zp_all(3,:)>1.28)),resp_ind));
    ind_comp(4) = length(intersect(intersect(find(Zc_all(4,:)>1.28),find(Zc_all(4,:)-Zp_all(4,:)>1.28)),resp_ind));

    avg_classification(1) = mean(ind_patt(:))/length(resp_ind);
    avg_classification(2) = mean(ind_comp(:))/length(resp_ind);

    ind_patt1 = intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28)),resp_ind);
    ind_patt2 = intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),resp_ind);
    ind_patt3 = intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(3,:)>1.28)),resp_ind);
    ind_patt4 = intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),resp_ind);
    ind_comp1 = intersect(intersect(find(Zc_all(1,:)>1.28),find(Zc_all(1,:)-Zp_all(1,:)>1.28)),resp_ind);
    ind_comp2 = intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),resp_ind);
    ind_comp3 = intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(3,:)-Zp_all(3,:)>1.28)),resp_ind);
    ind_comp4 = intersect(intersect(find(Zc_all(4,:)>1.28),find(Zc_all(4,:)-Zp_all(4,:)>1.28)),resp_ind);
    
    max_classification(1) = length(unique([ind_patt1;ind_patt2;ind_patt3;ind_patt4]))/length(resp_ind);
    max_classification(2) = length(unique([ind_comp1;ind_comp2;ind_comp3;ind_comp4]))/length(resp_ind);


    pie0 = [ind_patt(1)/length(resp_ind) ind_comp(1)/length(resp_ind) (1-(ind_patt(1)/length(resp_ind))-(ind_comp(1)/length(resp_ind)))];
    pie90 = [ind_patt(2)/length(resp_ind) ind_comp(2)/length(resp_ind) (1-(ind_patt(2)/length(resp_ind))-(ind_comp(2)/length(resp_ind)))];
    pie270 = [ind_patt(4)/length(resp_ind) ind_comp(4)/length(resp_ind) (1-(ind_patt(4)/length(resp_ind))-(ind_comp(4)/length(resp_ind)))];

    ind_patt_cells0 = intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28)),resp_ind);
    ind_patt_cells0_pattAt90 = intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_compAt90 = intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_pattAt180 = intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(3,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_compAt180 = intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(3,:)-Zp_all(3,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_pattAt270 = intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_compAt270 = intersect(intersect(find(Zc_all(4,:)>1.28),find(Zc_all(4,:)-Zp_all(4,:)>1.28)),ind_patt_cells0);
    

    ind_patt_cells270 = intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),resp_ind);
    ind_patt_cells270_pattAt90 = intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_compAt90 = intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_pattAt180 = intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(3,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_compAt180 = intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(3,:)-Zp_all(3,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_pattAt0 = intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(4,:)-Zc_all(1,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_compAt0 = intersect(intersect(find(Zc_all(1,:)>1.28),find(Zc_all(4,:)-Zp_all(1,:)>1.28)),ind_patt_cells270);

    pie_patt270at90 = [length(ind_patt_cells270_pattAt90)/length(ind_patt_cells270) length(ind_patt_cells270_compAt90)/length(ind_patt_cells270) (1-(length(ind_patt_cells270_pattAt90)/length(ind_patt_cells270))-( length(ind_patt_cells270_compAt90)/length(ind_patt_cells270)))];
    pie_patt270at180 = [length(ind_patt_cells270_pattAt180)/length(ind_patt_cells270) length(ind_patt_cells270_compAt180)/length(ind_patt_cells270) (1-(length(ind_patt_cells270_pattAt180)/length(ind_patt_cells270))-( length(ind_patt_cells270_compAt180)/length(ind_patt_cells270)))];
    pie_patt270at0 = [length(ind_patt_cells270_pattAt0)/length(ind_patt_cells270) length(ind_patt_cells270_compAt0)/length(ind_patt_cells270) (1-(length(ind_patt_cells270_pattAt0)/length(ind_patt_cells270))-( length(ind_patt_cells270_compAt0)/length(ind_patt_cells270)))];

    pie_patt0atOthers(1,:) = [pie_patt270at90];
    pie_patt0atOthers(2,:) = [pie_patt270at180];
    pie_patt0atOthers(3,:) = [pie_patt270at0];
    pie_patt0atOtherAvg = mean(pie_patt0atOthers,1);
    
    if iA == 1
        figure(2);
        subplot 331
            piechart(pie270)
            %subtitle('pie 0')
        subplot 332
            piechart(pie_patt0atOtherAvg)
            %subtitle('patt 0 at others avg')
        subplot 333
            piechart([max_classification(1) 1-max_classification(1)])
            %subtitle('max')
            movegui('center')
        sgtitle([leg_str{iA}])
    else
        figure(3);
        subplot 331
            piechart(pie270)
            %subtitle('pie 0')
        subplot 332
            piechart(pie_patt0atOtherAvg)
            %subtitle('patt 0 at others avg')
        subplot 333
            piechart([max_classification(1) 1-max_classification(1)])
            %subtitle('max')
            movegui('center')
        sgtitle([leg_str{iA}])
    end

end

%% add marm data


clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = ([base '\Analysis\Neuropixel\CrossOri\randDirFourPhase\summaries']);
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandDirRandPhaseSummary');
svName = 'randDirFourPhase_CrossOri';

load(fullfile(summaryDir,[svName '_Summary_V1_MAR.mat']))
resp_ind = intersect(intersect(sig_stim,sig_dir),find(DSI_all>0.5));


  nCellsResp = length(resp_ind);
    neuron_classes=[];
    clear ind_patt ind_comp

    for ip = 1:4
        ind_patt{ip} = intersect(find(Zp_all(ip,:)>1.28),find(Zp_all(ip,:)-Zc_all(ip,:)>1.28));
        ind_comp{ip} = intersect(find(Zc_all(ip,:)>1.28),find(Zc_all(ip,:)-Zp_all(ip,:)>1.28));
        ind_uncl{ip} = setdiff(1:length(Zc_all), [ind_patt{ip} ind_comp{ip}]);

        % Assign label for each neuron
        for ic = 1:nCellsResp
            i = resp_ind(ic);
            if ismember(i, ind_patt{ip})
                neuron_classes{ic}(ip) = 'p';
            elseif ismember(i, ind_comp{ip})
                neuron_classes{ic}(ip) = 'c';
            elseif ismember(i, ind_uncl{ip})
                neuron_classes{ic}(ip) = 'r';
            else
                neuron_classes{ic}(ip) = 'x'; % fallback in case itcs none
                error(['cell ' num2str(i) ' unclassed at phase=' num2str(ip)])
            end
        end
    end
   
    % Make lists for categories
    truly_p = [];
    truly_c = [];
    truly_r = [];
    switch_pc = [];
    switch_pr = [];
    switch_cr = [];
    switch_pcr = [];

    for i = 1:nCellsResp
        classes = unique(neuron_classes{i});
        classes(classes == 'x') = []; % remove unknowns if any
    
        switch length(classes)
            case 1
                switch classes
                    case 'p'
                        truly_p(end+1) = i;
                    case 'c'
                        truly_c(end+1) = i;
                    case 'r'
                        truly_r(end+1) = i;
                end
            case 2
                if all(ismember(classes, ['p' 'c']))
                    switch_pc(end+1) = i;
                elseif all(ismember(classes, ['p' 'r']))
                    switch_pr(end+1) = i;
                elseif all(ismember(classes, ['c' 'r']))
                    switch_cr(end+1) = i;
                end
            case 3
                switch_pcr(end+1) = i;
        end
    end

   % Labels and counts
    categories = {'Truly P', 'Truly C', 'Truly R', 'Switch P/C', 'Switch P/R', 'Switch C/R', 'Switch P/C/R'};
    counts = [length(truly_p), length(truly_c), length(truly_r), [length(switch_pc)+ length(switch_pr)+length(switch_cr)+length(switch_pcr)]];
    countsAll = [length(truly_p), length(truly_c), length(truly_r), ...
              length(switch_pc), length(switch_pr), length(switch_cr), length(switch_pcr)];
    
    % Convert to percentages
    percentages = 100 * counts / nCellsResp;
    
    % Plot
    figure(1);
    subplot(2,3,3)
    bar(percentages);
    set(gca, 'XTickLabel', categories, 'XTickLabelRotation', 45);
    ylabel('percentage of cells');
    ylim([0 100]);
    subtitle(['marm V1 n=' num2str(nCellsResp)])
    movegui('center');
    hold on
    sgtitle('Cell classification switching: p=pattern, c=component, r=unclassified');


    % pie charts
    clear ind_patt ind_comp

    ind_patt(1) = length(intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28)),resp_ind));
    ind_patt(2) = length(intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),resp_ind));
    ind_patt(3) = length(intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(2,:)>1.28)),resp_ind));
    ind_patt(4) = length(intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),resp_ind));
    ind_comp(1) = length(intersect(intersect(find(Zc_all(1,:)>1.28),find(Zc_all(1,:)-Zp_all(1,:)>1.28)),resp_ind));
    ind_comp(2) = length(intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),resp_ind));
    ind_comp(3) = length(intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(3,:)-Zp_all(3,:)>1.28)),resp_ind));
    ind_comp(4) = length(intersect(intersect(find(Zc_all(4,:)>1.28),find(Zc_all(4,:)-Zp_all(4,:)>1.28)),resp_ind));

    avg_classification(1) = mean(ind_patt(:))/length(resp_ind);
    avg_classification(2) = mean(ind_comp(:))/length(resp_ind);

    ind_patt1 = intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28)),resp_ind);
    ind_patt2 = intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),resp_ind);
    ind_patt3 = intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(3,:)>1.28)),resp_ind);
    ind_patt4 = intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),resp_ind);
    ind_comp1 = intersect(intersect(find(Zc_all(1,:)>1.28),find(Zc_all(1,:)-Zp_all(1,:)>1.28)),resp_ind);
    ind_comp2 = intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),resp_ind);
    ind_comp3 = intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(3,:)-Zp_all(3,:)>1.28)),resp_ind);
    ind_comp4 = intersect(intersect(find(Zc_all(4,:)>1.28),find(Zc_all(4,:)-Zp_all(4,:)>1.28)),resp_ind);
    
    max_classification(1) = length(unique([ind_patt1;ind_patt2;ind_patt3;ind_patt4]))/length(resp_ind);
    max_classification(2) = length(unique([ind_comp1;ind_comp2;ind_comp3;ind_comp4]))/length(resp_ind);

    pie0 = [ind_patt(1)/length(resp_ind) ind_comp(1)/length(resp_ind) (1-(ind_patt(1)/length(resp_ind))-(ind_comp(1)/length(resp_ind)))];
    pie90 = [ind_patt(2)/length(resp_ind) ind_comp(2)/length(resp_ind) (1-(ind_patt(2)/length(resp_ind))-(ind_comp(2)/length(resp_ind)))];
    pie270 = [ind_patt(4)/length(resp_ind) ind_comp(4)/length(resp_ind) (1-(ind_patt(4)/length(resp_ind))-(ind_comp(4)/length(resp_ind)))];

    ind_patt_cells0 = intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(1,:)-Zc_all(1,:)>1.28)),resp_ind);
    ind_patt_cells0_pattAt90 = intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_compAt90 = intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_pattAt180 = intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(3,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_compAt180 = intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(3,:)-Zp_all(3,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_pattAt270 = intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),ind_patt_cells0);
    ind_patt_cells0_compAt270 = intersect(intersect(find(Zc_all(4,:)>1.28),find(Zc_all(4,:)-Zp_all(4,:)>1.28)),ind_patt_cells0);

    ind_patt_cells270 = intersect(intersect(find(Zp_all(4,:)>1.28),find(Zp_all(4,:)-Zc_all(4,:)>1.28)),resp_ind);
    ind_patt_cells270_pattAt90 = intersect(intersect(find(Zp_all(2,:)>1.28),find(Zp_all(2,:)-Zc_all(2,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_compAt90 = intersect(intersect(find(Zc_all(2,:)>1.28),find(Zc_all(2,:)-Zp_all(2,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_pattAt180 = intersect(intersect(find(Zp_all(3,:)>1.28),find(Zp_all(3,:)-Zc_all(3,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_compAt180 = intersect(intersect(find(Zc_all(3,:)>1.28),find(Zc_all(3,:)-Zp_all(3,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_pattAt0 = intersect(intersect(find(Zp_all(1,:)>1.28),find(Zp_all(4,:)-Zc_all(1,:)>1.28)),ind_patt_cells270);
    ind_patt_cells270_compAt0 = intersect(intersect(find(Zc_all(1,:)>1.28),find(Zc_all(4,:)-Zp_all(1,:)>1.28)),ind_patt_cells270);

    pie_patt270at90 = [length(ind_patt_cells270_pattAt90)/length(ind_patt_cells270) length(ind_patt_cells270_compAt90)/length(ind_patt_cells270) (1-(length(ind_patt_cells270_pattAt90)/length(ind_patt_cells270))-( length(ind_patt_cells270_compAt90)/length(ind_patt_cells270)))];
    pie_patt270at180 = [length(ind_patt_cells270_pattAt180)/length(ind_patt_cells270) length(ind_patt_cells270_compAt180)/length(ind_patt_cells270) (1-(length(ind_patt_cells270_pattAt180)/length(ind_patt_cells270))-( length(ind_patt_cells270_compAt180)/length(ind_patt_cells270)))];
    pie_patt270at0 = [length(ind_patt_cells270_pattAt0)/length(ind_patt_cells270) length(ind_patt_cells270_compAt0)/length(ind_patt_cells270) (1-(length(ind_patt_cells270_pattAt0)/length(ind_patt_cells270))-( length(ind_patt_cells270_compAt0)/length(ind_patt_cells270)))];

    pie_patt0atOthers(1,:) = [pie_patt270at90];
    pie_patt0atOthers(2,:) = [pie_patt270at180];
    pie_patt0atOthers(3,:) = [pie_patt270at0];
    pie_patt0atOtherAvg = mean(pie_patt0atOthers,1);
    
        figure(4);
        subplot 331
            piechart(pie270)
            %subtitle('pie 0')
        subplot 332
            piechart(pie_patt0atOtherAvg)
            %subtitle('patt 0 at others avg')
        subplot 333
            piechart([max_classification(1) 1-max_classification(1)])
            %subtitle('max')
            movegui('center')
        sgtitle('marmoset')

%%
figure(1); 
    print(fullfile(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara\Figures\CrossOri_meetings\'], ['TrulyPattCompUncl_BarCharts.pdf']), '-dpdf','-fillpage')
figure(2); 
    print(fullfile(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara\Figures\CrossOri_meetings\'], ['V1_mouse_layer4_PieCharts.pdf']), '-dpdf','-fillpage')
figure(3); 
    print(fullfile(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara\Figures\CrossOri_meetings\'], ['V1_mouse_layer23_PieCharts.pdf']), '-dpdf','-fillpage')
figure(4); 
    print(fullfile(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara\Figures\CrossOri_meetings\'], ['V1_marmoset_PieCharts.pdf']), '-dpdf','-fillpage')

        