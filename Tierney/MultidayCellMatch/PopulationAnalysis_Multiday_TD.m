<<<<<<< HEAD
%% Experiment info
clear all; clear global; close all
clc

%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
fn_analysis = fullfile(fn_base, 'home\Tierney\Analysis\2P');
fn_multi = fullfile(fn_base, 'home\Tierney\Analysis\2P\MultidayAnalysis');
fn_pop = fullfile(fn_base, 'home\Tierney\Analysis\2P\PopulationAnalysis');

ds = 'ExperimentData_TD'; % dataset info 
eval(ds)

sess_list = [2, 4, 8, 11];% enter all the sessions you want to pool
nSess = length(sess_list);
 
% sess_title = string(sess_list(1));
% for iSess = 2:nSess
%     sess_title = strcat(sess_title,'_',string(sess_list(iSess)));
% end

for iSess = 1:nSess
    % Baseline (1) and post-MD (2) days
    day_id(2) = sess_list(iSess);
    day_id(1) = expt(day_id(2)).matchday_baseline;

    % % Baseline (1), post-MD (2), and recovery (3) days
    % day_id(2) = 2;
    % day_id(1) = expt(day_id(2)).matchday_baseline;
    % day_id(3) = expt(day_id(2)).matchday_recovery;

    nd = length(day_id);

    %Specific experiment information
    mouse = expt(day_id(1)).mouse;

    expt(day_id(2)).multiday_time_days>0
    time_str = [num2str(expt(day_id(2)).multiday_time_days) 'Days'];
    fn_match = ['multiday_' time_str,'_',expt(day_id(2)).experiment];

    for id = 1:nd
        date = expt(day_id(id)).date;
        runs = expt(day_id(id)).stimruns;
        nrun = length(runs);
        contra = strcmp(expt(day_id(id)).eye_str,'Contra'); % 1 is contra eye open; 0 is ipsi eye open

        run_str = catRunName(runs, nrun);
        eye_str = expt(day_id(id)).eye_str;

        datemouse = [date '_' mouse];
        datemouserun = [date '_' mouse '_' run_str];

%         % Load respData and ODI
%         fn_stimData = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun '_stimData_Ori.mat']);
% %         fn_respData = fullfile(fn_multi, mouse, fn_match, [mouse, '_respData_multiday.mat']);
% %         fn_ODI = fullfile(fn_multi, mouse, fn_match, [mouse, '_ODI_multiday.mat']);
% %         fn_oriResp = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun, '_oriResp' '.mat']);
% 
%         temp_stimData = load(fn_stimData{iSess});
% %         temp_respData = load(fn_respData);
% %         temp_ODI = load(fn_ODI);
% %         temp_oriResp = load(fn_oriResp);
% 
%         stimData{id} = temp_stimData;
% %         respData{id} = temp_respData;
% %         ODI{id} = temp_ODI;
% %         oriResp{id} = temp_oriResp;
% 
%         clear temp_stimData
% %         clear temp_respData
% %         clear temp_ODI
% %         clear temp_oriResp
    end

    % Load stimData, timecourses, and input
    sess_TCs{iSess} = load(fullfile(fn_multi, mouse, fn_match, ['timecourses.mat']));
    sess_input{iSess} = load(fullfile(fn_multi, mouse, fn_match, ['input.mat']));
    sess_stimData{iSess} = load(fullfile(fn_multi, mouse, fn_match, [mouse '_stimData_multiday.mat']));
end

%% Extract variables of interest from the cell arrays loaded above

% Use this section if stimData values are different for each day; would need
% to add (id) in for loops after these variables.

% nOn = zeros(1,nd); 
% nOff = zeros(1,nd); 
% ntrials = zeros(1,nd); 
% Eyes = zeros(2,nd); 
% nEye = zeros(1,nd); 
% Dirs = cell(1,nd); 
% nDirs = zeros(1,nd); 

%  for id = 1:nd
%     nOn(id) = stimData{id}.nOn;
%     nOff(id) = stimData{id}.nOff;
%     ntrials(id) = stimData{id}.ntrials;
%     
%     tContra{id} = stimData{id}.tContra;
%     Eyes(:,id) = stimData{id}.Eyes;
%     nEye(id) = stimData{id}.nEye;
% 
%     tDir{id} = stimData{id}.tDir;
%     Dirs{id} = stimData{id}.Dirs;
%     nDirs(id) = stimData{id}.nDirs;
%  end

tContra = cell(1,nd);
tOri = cell(1,nd);

nOn = sess_stimData{1}.nOn;
nOff = sess_stimData{1}.nOff;

Eyes = sess_stimData{1}.Eyes;
nEye = sess_stimData{1}.nEye;

Oris = sess_stimData{1}.Oris;
nOris = sess_stimData{1}.nOris;

for iSess = 1:nSess
    for id = 1:nd
        ntrials(iSess,id) = sess_stimData{iSess}.ntrials(id);
        tContra{iSess, id} = sess_stimData{iSess}.tContra{id};
        tOri{iSess, id} = sess_stimData{iSess}.tOri{id};
        cellTCs_match{iSess, id} = sess_TCs{iSess}.cellTCs_match{id};
    end
end

%% looking at time courses: average across all trials
% 
% data_dfof_trial = cell(1,nd);
% 
% % In the below plots, circshift (-50) was used, so stimOn starts at 10
% % frames.
% figure
% for id = 1:nd
%     nCells = size(cellTCs_match{id},2);
% 
%     data_tc_trial{id} = reshape(cellTCs_match{id}, [nOn+nOff,ntrials(id),nCells]);
%     data_f_trial{id} = mean(data_tc_trial{id}(nOff/2:nOff,:,:),1);
%     data_dfof_trial{id} = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial{id}, data_f_trial{id}), data_f_trial{id});
% 
%     %looking at data with np subtracted
%     tc_cell_avrg{id} = mean(data_dfof_trial{id},3);%average per cells, one row per trial
%     tc_trial_avrg{id} = squeeze(mean(data_dfof_trial{id},2));%average over trials, one row per cell
%     tc_cell_trial_avrg{id} = mean(tc_cell_avrg{id},2);%average over trials and cells
% 
%     subplot(2,1,id)
%     plot(circshift(tc_trial_avrg{id},-50), 'LineWidth',.005);
%     hold on;
%     plot(circshift(tc_cell_trial_avrg{id},-50), 'LineWidth',2, 'color','k');
%     hold on;
%     title(['Timecourses: ' num2str(expt(id).experiment)]); % Timecourses with neuropil subtracted
%     xlabel('Frames')
%     ylabel('df/f')
%     hold off
% end
% 
% print(fullfile(fn_pop, ['Timecourses_population.pdf']),'-dpdf','-bestfit')

%% get tuning data
base_win = nOff/2:nOff;
resp_win = nOff+5:nOff+nOn;

for iSess = 1:nSess
    for id = 1:nd
        nCells(iSess,id) = size(cellTCs_match{iSess,id},2);
        data_trial{iSess,id} = reshape(cellTCs_match{iSess,id}, [nOn+nOff ntrials(iSess,id) nCells(iSess,id)]);
        data_f{iSess,id} = mean(data_trial{iSess,id}(base_win,:,:),1);
        data_dfof{iSess,id} = bsxfun(@rdivide,bsxfun(@minus,data_trial{iSess,id},data_f{iSess,id}),data_f{iSess,id});
    end
end

resp_cell_Ori = cell(nSess,nEye,nOris,nd);
base_cell_Ori = cell(nSess,nEye,nOris,nd);
data_dfof_Ori = cell(nSess,nSess,nd);
h_Ori = cell(nSess,nd);
p_Ori = cell(nSess,nd);

for iSess = 1:nSess
    for id = 1:nd
        data_dfof_Ori_temp = zeros(nCells(iSess,id),nOris,nEye,2);
        h_Ori_temp = zeros(nOris,nCells(iSess,id),nEye);
        p_Ori_temp = zeros(nOris,nCells(iSess,id),nEye);

        for iEye = 1:nEye
            ind_eye = find(tContra{id} == Eyes(iEye));
            for iOri = 1:nOris
                ind_Ori = find(tOri{iSess,id} == Oris(iOri));
                ind = intersect(ind_eye,ind_Ori);
                resp_cell_Ori{iSess,iEye,iOri,id} = squeeze(mean(data_dfof{iSess,id}(resp_win,ind,:),1));
                base_cell_Ori{iSess,iEye,iOri,id} = squeeze(mean(data_dfof{iSess,id}(base_win,ind,:),1));
                [h_Ori_temp(iOri,:,iEye), p_Ori_temp(iOri,:,iEye)] = ttest(resp_cell_Ori{iSess,iEye,iOri,id},base_cell_Ori{iSess,iEye,iOri,id},'tail','right','alpha',0.05./((nOris*nEye)-1));
                data_dfof_Ori_temp(:,iOri,iEye,1) = squeeze(mean(mean(data_dfof{iSess,id}(resp_win,ind,:),1),2));
                data_dfof_Ori_temp(:,iOri,iEye,2) = squeeze(std(mean(data_dfof{iSess,id}(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
            end
        end

        h_Ori{iSess,id} = h_Ori_temp;
        p_Ori{iSess,id} = p_Ori_temp;
        data_dfof_Ori{iSess,id} = data_dfof_Ori_temp;

        h_all_Ori{iSess,id} = squeeze(sum(h_Ori{iSess,id},1));
        resp_ind{iSess,id} = find(sum(h_all_Ori{iSess,id},2));
        ipsi_resp_ind{iSess,id} = find(h_all_Ori{iSess,id}(:,1));
        contra_resp_ind{iSess,id} = find(h_all_Ori{iSess,id}(:,2));
    end 
end

clear h_Ori_temp p_Ori_temp data_dfof_Ori_temp

save(fullfile(fn_pop, ['respData_population.mat']), 'h_Ori', 'resp_ind', 'ipsi_resp_ind', 'contra_resp_ind', 'resp_cell_Ori', 'base_cell_Ori', 'data_dfof_Ori','base_win','resp_win','data_dfof')

figure
color_choices = [1,0,0; 0,1,0; 0,0,1; 0,1,1];
sz_dots = 15
for iEye = 1:nEye
    subplot(2,1,iEye)
    for iSess = 1:nSess
    %     C = linspace(1,10,length(h_all_Ori{2}(:,iEye)));
    %     swarmchart((h_all_Ori{1}(:,iEye)),(h_all_Ori{2}(:,iEye)),[],C)
        swarmchart((h_all_Ori{iSess,1}(:,iEye)),(h_all_Ori{iSess,2}(:,iEye)), sz_dots,...
            'MarkerFaceColor', color_choices(iSess,:),...
            'MarkerEdgeColor', color_choices(iSess,:),...
            'MarkerFaceAlpha', 0.5,...
            'MarkerEdgeAlpha', 0.5)
        axis square
        title(['Significant orientations: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
        xlabel([num2str(expt(1).experiment)])
        ylabel([num2str(expt(2).experiment)])
        refline(1)
        hold on 
    end
end

print(fullfile(fn_pop, ['sigOrisScatter_population.pdf']),'-dpdf','-bestfit')

%% get ODI

for iSess = 1:nSess
    for id = 1:nd
        [contra_resp{iSess,id} max_ind_contra{iSess,id}] = max(data_dfof_Ori{iSess,id}(:,:,find(Eyes),1),[],2);
        max_Ori_contra{iSess,id} = Oris(max_ind_contra{iSess,id});
        [ipsi_resp{iSess,id} max_ind_ipsi{iSess,id}] = max(data_dfof_Ori{iSess,id}(:,:,find(~Eyes),1),[],2); 
        max_Ori_ipsi{iSess,id} = Oris(max_ind_ipsi{iSess,id});

        real_contra_resp{iSess,id} = contra_resp{iSess,id};
        real_contra_ind{iSess,id} = find(real_contra_resp{iSess,id} < 0);
        real_contra_resp{iSess,id}(real_contra_ind{iSess,id}) = 0;

        real_ipsi_resp{iSess,id} = ipsi_resp{iSess,id};
        real_ipsi_ind{iSess,id} = find(real_ipsi_resp{iSess,id} < 0);
        real_ipsi_resp{iSess,id}(real_ipsi_ind{iSess,id}) = 0;

        ODI{iSess,id} = (real_contra_resp{iSess,id}-real_ipsi_resp{iSess,id})./(real_contra_resp{iSess,id}+real_ipsi_resp{iSess,id});
    end
end

%scatter of max response to contra and ipsi
figure
for id = 1:nd
    for iSess = 1:nSess
    subplot(2,2,id)
    contra_resp_any{iSess,id} = max(data_dfof_Ori{iSess,id}(resp_ind{iSess,id},:,find(Eyes),1),[],2);
    ipsi_resp_any{iSess,id} = max(data_dfof_Ori{iSess,id}(resp_ind{iSess,id},:,find(~Eyes),1),[],2);
%     C = linspace(1,10,length(contra_resp_any{iSess,id}(:,1)));
%     scatter(contra_resp_any{iSess,id},ipsi_resp_any{iSess,id},[],C)
    swarmchart(contra_resp_any{iSess,id},ipsi_resp_any{iSess,id}, sz_dots,...
            'MarkerFaceColor', color_choices(iSess,:),...
            'MarkerEdgeColor', color_choices(iSess,:),...
            'MarkerFaceAlpha', 0.5,...
            'MarkerEdgeAlpha', 0.5)
    hold on
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    axis square
    xlabel('Max df/f (Contra)')
    ylabel('Max df/f (Ipsi)')
    xlim([0.001 1])
    ylim([0.001 1])
    refline(1)
    title(['Responsive cells: ' num2str(expt(id).experiment)])
    end
end

% All ODI

for iSess = 1:nSess
    subplot(2,2,id+1)
    swarmchart(ODI{iSess,1}, ODI{iSess,2}, sz_dots,...
            'MarkerFaceColor', color_choices(iSess,:),...
            'MarkerEdgeColor', color_choices(iSess,:),...
            'MarkerFaceAlpha', 0.5,...
            'MarkerEdgeAlpha', 0.5)
    hold on
    axis square
    xlim([-1 1])
    ylim([-1 1])
    refline(1)
    xlabel([num2str(expt(1).experiment)])
    ylabel([num2str(expt(2).experiment)])
    title('ODI')
end

% Mean ODI
mean_ODI = zeros(nSess,nd);

for iSess = 1:nSess 
    for id = 1:nd
        mean_ODI(iSess,id) = nanmean(ODI{iSess,id});
    end
end

for iSess = 1:nSess
    subplot(2,2,id+2)
    swarmchart(mean_ODI(iSess,1), mean_ODI(iSess,2), sz_dots,...
            'MarkerFaceColor', color_choices(iSess,:),...
            'MarkerEdgeColor', color_choices(iSess,:),...
            'MarkerFaceAlpha', 0.5,...
            'MarkerEdgeAlpha', 0.5)
    hold on
    axis square
    xlim([-1 1])
    ylim([-1 1])
    refline(1)
    xlabel([num2str(expt(1).experiment)])
    ylabel([num2str(expt(2).experiment)])
    title('Mean ODI')
end

print(fullfile(fn_pop, ['OD_population.pdf']),'-dpdf','-bestfit')

save(fullfile(fn_pop, ['ODI_population.mat']), 'ODI', 'real_contra_resp', 'real_ipsi_resp', 'max_Ori_contra', 'max_Ori_ipsi')

%% Population averages for ODI

% % Do this for recovery, but histogram
% swarmchart(day_id(1:3),mean_ODI)
% axis square
% xlabel('Timepoint')
% set(gca,'XTick',3:5,'XTickLabel',{'Pre','Post','Rec'})
% ylabel('Mean ODI')
% ylim([-1 1])
% title('Mean ODI across days')

%% Population average for contra/ipsi inputs

% Mean contra/ipsi inputs
mean_ipsi = zeros(nSess,nd);
mean_contra = zeros(nSess,nd);

for iSess = 1:nSess
    for id = 1:nd
         mean_ipsi(iSess,id) = mean(real_ipsi_resp{iSess,id});
         mean_contra(iSess,id) = mean(real_contra_resp{iSess,id});
         mean_resp(iSess,id,:) = [mean_ipsi(iSess,id); mean_contra(iSess,id)]
    end
end

figure
for iSess = 1:nSess
    scatter([1 2 3 4],reshape(mean_resp(iSess,:,:),[1,4]), sz_dots,...
            'MarkerFaceColor', color_choices(iSess,:),...
            'MarkerEdgeColor', color_choices(iSess,:),...
            'MarkerFaceAlpha', 0.5,...
            'MarkerEdgeAlpha', 0.5)
    hold on
    axis square
    xlabel('Timepoint')
    set(gca,'XTick',1:4,'XTickLabel',{'Pre: Ipsi', 'Post: Ipsi','Pre: Contra','Post: Contra'})
    ylabel('Mean response (df/f)')
    title(['Ipsi v. contra response'])
    axis([0.5 4.5 0 0.35])
end

print(fullfile(fn_pop, ['EyeRespHist_population.pdf']),'-dpdf','-bestfit')

save(fullfile(fn_pop, ['EyeResp_population.mat']), 'mean_ipsi', 'mean_contra')

%% Need to finish editing this code
% %% von mises
% 
%     b_ori = cell(1,nd);
%     k1_ori = cell(1,nd);
%     R1_ori = cell(1,nd);
%     u1_ori = cell(1,nd);
%     R_square_ori = cell(1,nd);
%     sse_ori = cell(1,nd);
%     stim_OSI = cell(1,nd);
%     theta_hires = deg2rad(0:180);
%     y_fit = cell(1,nd);
% 
% for id = 1:nd    
%     for iEye = 1:nEye
%         for iCell = 1:nCells(iSess,id)
%             data = [data_dfof_Ori{id}(iCell,:,iEye) data_dfof_Ori{id}(iCell,1,iEye)];
%             theta = [deg2rad(Oris) pi];
%             [b_ori{id}(iEye,iCell),k1_ori{id}(iEye,iCell),R1_ori{id}(iEye,iCell),u1_ori{id}(iEye,iCell),sse_ori{id}(iEye,iCell),R_square_ori{id}(iEye,iCell)] ...
%                 = miaovonmisesfit_ori(theta,data);
%             [max_val{id} max_ind{id}] = max(data_dfof_Ori{id}(iCell,:,iEye),[],2);
%             null_ind{id} = max_ind{id}+(nOris./2);
%             null_ind{id}(find(null_ind{id}>nOris)) = null_ind{id}(find(null_ind{id}>nOris))-nOris;
%             min_val{id} = data_dfof_Ori{id}(iCell,null_ind{id},iEye);
%             if min_val{id}<0
%                 min_val{id} = 0;
%             end
%             stim_OSI{id}(1,iCell) = (max_val{id}-min_val{id})./(max_val{id}+min_val{id});
%             y_fit{id}(:,iEye,iCell) = b_ori{id}(iEye,iCell) + R1_ori{id}(iEye,iCell) .* exp(k1_ori{id}(iEye,iCell).*(cos(2.*(theta_hires-u1_ori{id}(iEye,iCell)))-1));
%         end
%     end
%     
%     [yfit_max{id}, yfit_max_ind{id}] = max(y_fit{id}(:,:,:),[],1);
%     prefOri_yfit{id} = squeeze(theta_hires(yfit_max_ind{id}));
% end   
% save(fullfile(fn_pop, ['oriResp_population.mat']), 'data_dfof_Ori','base_win','resp_win','h_Ori', 'b_ori', 'k1_ori', 'R1_ori', 'u1_ori', 'R_square_ori', 'sse_ori','stim_OSI', 'yfit_max', 'yfit_max_ind', 'prefOri_yfit') 
% 
% %% Population average for tuning.
% 
% % Plot tuning widths in cumulative distribution plot
% figure
% for iEye = 1:nEye
%     subplot(2,2,iEye)
%     for id = 1:nd
%         cdfplot(k1_ori{id}(iEye,:));
%         hold on
%     end
%     axis square
%     xlabel('Tuning width')
%     ylabel('Fraction of cells')
%     title(['Tuning width: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
% end
% legend('Pre','Post','Rec','Location','Southeast');
% 
% % Compare tuning widths of individual cells -- plot in scatterplots
% % (baseline v. MD)
% for iEye = 1:nEye
%     subplot(2,2,iEye+2)
%     scatter(k1_ori{1}(iEye,:),k1_ori{2}(iEye,:));
%     axis square
%     title(['Tuning width: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
%     xlabel([num2str(expt(1).experiment)])
%     ylabel([num2str(expt(2).experiment)])
%     refline(1)
%     hold on   
% end
% 
% print(fullfile(fn_pop, ['TuningPlots_population.pdf']),'-dpdf','-bestfit')
% 
% %% Compare dfof between pre- and post-MD
% 
% % EDIT: FINISH THIS SECTION
% figure
% for iEye = nEye
%     
% end
% 
% temp_data = data_dfof{1};
% temp_max = max(reshape(temp_data,[],size(temp_data,3)),[],1);
% 
% % EDIT: PUT TEMP_MAX IN HERE INSTEAD????
% [contra_resp{id} max_ind_contra{id}] = max(data_dfof_Ori{id}(:,:,find(Eyes),1),[],2);
% 
% subplot(2,1,iEye)
% scatter((data_dfof_max{1}(:,iEye)),(data_dfof_max{2}(:,iEye)))
% axis square
% title(['Max dfof: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
% xlabel([num2str(expt(1).experiment)])
% ylabel([num2str(expt(2).experiment)])
% refline(1)
=======
%% Experiment info
clear all; clear global; close all
clc

%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
fn_analysis = fullfile(fn_base, 'home\Tierney\Analysis\2P');
fn_multi = fullfile(fn_base, 'home\Tierney\Analysis\2P\MultidayAnalysis');
fn_pop = fullfile(fn_base, 'home\Tierney\Analysis\2P\PopulationAnalysis');

ds = 'ExperimentData_TD'; % dataset info 
eval(ds)

sess_list = [2, 4, 8, 11];% enter all the sessions you want to pool
nSess = length(sess_list);
 
% sess_title = string(sess_list(1));
% for iSess = 2:nSess
%     sess_title = strcat(sess_title,'_',string(sess_list(iSess)));
% end

for iSess = 1:nSess
    % Baseline (1) and post-MD (2) days
    day_id(2) = sess_list(iSess);
    day_id(1) = expt(day_id(2)).matchday_baseline;

    % % Baseline (1), post-MD (2), and recovery (3) days
    % day_id(2) = 2;
    % day_id(1) = expt(day_id(2)).matchday_baseline;
    % day_id(3) = expt(day_id(2)).matchday_recovery;

    nd = length(day_id);

    %Specific experiment information
    mouse = expt(day_id(1)).mouse;

    expt(day_id(2)).multiday_time_days>0
    time_str = [num2str(expt(day_id(2)).multiday_time_days) 'Days'];
    fn_match = ['multiday_' time_str,'_',expt(day_id(2)).experiment];

    for id = 1:nd
        date = expt(day_id(id)).date;
        runs = expt(day_id(id)).stimruns;
        nrun = length(runs);
        contra = strcmp(expt(day_id(id)).eye_str,'Contra'); % 1 is contra eye open; 0 is ipsi eye open

        run_str = catRunName(runs, nrun);
        eye_str = expt(day_id(id)).eye_str;

        datemouse = [date '_' mouse];
        datemouserun = [date '_' mouse '_' run_str];

%         % Load respData and ODI
%         fn_stimData = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun '_stimData_Ori.mat']);
% %         fn_respData = fullfile(fn_multi, mouse, fn_match, [mouse, '_respData_multiday.mat']);
% %         fn_ODI = fullfile(fn_multi, mouse, fn_match, [mouse, '_ODI_multiday.mat']);
% %         fn_oriResp = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun, '_oriResp' '.mat']);
% 
%         temp_stimData = load(fn_stimData{iSess});
% %         temp_respData = load(fn_respData);
% %         temp_ODI = load(fn_ODI);
% %         temp_oriResp = load(fn_oriResp);
% 
%         stimData{id} = temp_stimData;
% %         respData{id} = temp_respData;
% %         ODI{id} = temp_ODI;
% %         oriResp{id} = temp_oriResp;
% 
%         clear temp_stimData
% %         clear temp_respData
% %         clear temp_ODI
% %         clear temp_oriResp
    end

    % Load stimData, timecourses, and input
    sess_TCs{iSess} = load(fullfile(fn_multi, mouse, fn_match, ['timecourses.mat']));
    sess_input{iSess} = load(fullfile(fn_multi, mouse, fn_match, ['input.mat']));
    sess_stimData{iSess} = load(fullfile(fn_multi, mouse, fn_match, [mouse '_stimData_multiday.mat']));
end

%% Extract variables of interest from the cell arrays loaded above

% Use this section if stimData values are different for each day; would need
% to add (id) in for loops after these variables.

% nOn = zeros(1,nd); 
% nOff = zeros(1,nd); 
% ntrials = zeros(1,nd); 
% Eyes = zeros(2,nd); 
% nEye = zeros(1,nd); 
% Dirs = cell(1,nd); 
% nDirs = zeros(1,nd); 

%  for id = 1:nd
%     nOn(id) = stimData{id}.nOn;
%     nOff(id) = stimData{id}.nOff;
%     ntrials(id) = stimData{id}.ntrials;
%     
%     tContra{id} = stimData{id}.tContra;
%     Eyes(:,id) = stimData{id}.Eyes;
%     nEye(id) = stimData{id}.nEye;
% 
%     tDir{id} = stimData{id}.tDir;
%     Dirs{id} = stimData{id}.Dirs;
%     nDirs(id) = stimData{id}.nDirs;
%  end

tContra = cell(1,nd);
tOri = cell(1,nd);

nOn = sess_stimData{1}.nOn;
nOff = sess_stimData{1}.nOff;

Eyes = sess_stimData{1}.Eyes;
nEye = sess_stimData{1}.nEye;

Oris = sess_stimData{1}.Oris;
nOris = sess_stimData{1}.nOris;

for iSess = 1:nSess
    for id = 1:nd
        ntrials(iSess,id) = sess_stimData{iSess}.ntrials(id);
        tContra{iSess, id} = sess_stimData{iSess}.tContra{id};
        tOri{iSess, id} = sess_stimData{iSess}.tOri{id};
        cellTCs_match{iSess, id} = sess_TCs{iSess}.cellTCs_match{id};
    end
end

%% looking at time courses: average across all trials
% 
% data_dfof_trial = cell(1,nd);
% 
% % In the below plots, circshift (-50) was used, so stimOn starts at 10
% % frames.
% figure
% for id = 1:nd
%     nCells = size(cellTCs_match{id},2);
% 
%     data_tc_trial{id} = reshape(cellTCs_match{id}, [nOn+nOff,ntrials(id),nCells]);
%     data_f_trial{id} = mean(data_tc_trial{id}(nOff/2:nOff,:,:),1);
%     data_dfof_trial{id} = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial{id}, data_f_trial{id}), data_f_trial{id});
% 
%     %looking at data with np subtracted
%     tc_cell_avrg{id} = mean(data_dfof_trial{id},3);%average per cells, one row per trial
%     tc_trial_avrg{id} = squeeze(mean(data_dfof_trial{id},2));%average over trials, one row per cell
%     tc_cell_trial_avrg{id} = mean(tc_cell_avrg{id},2);%average over trials and cells
% 
%     subplot(2,1,id)
%     plot(circshift(tc_trial_avrg{id},-50), 'LineWidth',.005);
%     hold on;
%     plot(circshift(tc_cell_trial_avrg{id},-50), 'LineWidth',2, 'color','k');
%     hold on;
%     title(['Timecourses: ' num2str(expt(id).experiment)]); % Timecourses with neuropil subtracted
%     xlabel('Frames')
%     ylabel('df/f')
%     hold off
% end
% 
% print(fullfile(fn_pop, ['Timecourses_population.pdf']),'-dpdf','-bestfit')

%% get tuning data
base_win = nOff/2:nOff;
resp_win = nOff+5:nOff+nOn;

for iSess = 1:nSess
    for id = 1:nd
        nCells(iSess,id) = size(cellTCs_match{iSess,id},2);
        data_trial{iSess,id} = reshape(cellTCs_match{iSess,id}, [nOn+nOff ntrials(iSess,id) nCells(iSess,id)]);
        data_f{iSess,id} = mean(data_trial{iSess,id}(base_win,:,:),1);
        data_dfof{iSess,id} = bsxfun(@rdivide,bsxfun(@minus,data_trial{iSess,id},data_f{iSess,id}),data_f{iSess,id});
    end
end

resp_cell_Ori = cell(nSess,nEye,nOris,nd);
base_cell_Ori = cell(nSess,nEye,nOris,nd);
data_dfof_Ori = cell(nSess,nSess,nd);
h_Ori = cell(nSess,nd);
p_Ori = cell(nSess,nd);

for iSess = 1:nSess
    for id = 1:nd
        data_dfof_Ori_temp = zeros(nCells(iSess,id),nOris,nEye,2);
        h_Ori_temp = zeros(nOris,nCells(iSess,id),nEye);
        p_Ori_temp = zeros(nOris,nCells(iSess,id),nEye);

        for iEye = 1:nEye
            ind_eye = find(tContra{id} == Eyes(iEye));
            for iOri = 1:nOris
                ind_Ori = find(tOri{iSess,id} == Oris(iOri));
                ind = intersect(ind_eye,ind_Ori);
                resp_cell_Ori{iSess,iEye,iOri,id} = squeeze(mean(data_dfof{iSess,id}(resp_win,ind,:),1));
                base_cell_Ori{iSess,iEye,iOri,id} = squeeze(mean(data_dfof{iSess,id}(base_win,ind,:),1));
                [h_Ori_temp(iOri,:,iEye), p_Ori_temp(iOri,:,iEye)] = ttest(resp_cell_Ori{iSess,iEye,iOri,id},base_cell_Ori{iSess,iEye,iOri,id},'tail','right','alpha',0.05./((nOris*nEye)-1));
                data_dfof_Ori_temp(:,iOri,iEye,1) = squeeze(mean(mean(data_dfof{iSess,id}(resp_win,ind,:),1),2));
                data_dfof_Ori_temp(:,iOri,iEye,2) = squeeze(std(mean(data_dfof{iSess,id}(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
            end
        end

        h_Ori{iSess,id} = h_Ori_temp;
        p_Ori{iSess,id} = p_Ori_temp;
        data_dfof_Ori{iSess,id} = data_dfof_Ori_temp;

        h_all_Ori{iSess,id} = squeeze(sum(h_Ori{iSess,id},1));
        resp_ind{iSess,id} = find(sum(h_all_Ori{iSess,id},2));
        ipsi_resp_ind{iSess,id} = find(h_all_Ori{iSess,id}(:,1));
        contra_resp_ind{iSess,id} = find(h_all_Ori{iSess,id}(:,2));
    end 
end

clear h_Ori_temp p_Ori_temp data_dfof_Ori_temp

save(fullfile(fn_pop, ['respData_population.mat']), 'h_Ori', 'resp_ind', 'ipsi_resp_ind', 'contra_resp_ind', 'resp_cell_Ori', 'base_cell_Ori', 'data_dfof_Ori','base_win','resp_win','data_dfof')

figure
color_choices = [1,0,0; 0,1,0; 0,0,1; 0,1,1];
sz_dots = 15
for iEye = 1:nEye
    subplot(2,1,iEye)
    for iSess = 1:nSess
    %     C = linspace(1,10,length(h_all_Ori{2}(:,iEye)));
    %     swarmchart((h_all_Ori{1}(:,iEye)),(h_all_Ori{2}(:,iEye)),[],C)
        swarmchart((h_all_Ori{iSess,1}(:,iEye)),(h_all_Ori{iSess,2}(:,iEye)), sz_dots,...
            'MarkerFaceColor', color_choices(iSess,:),...
            'MarkerEdgeColor', color_choices(iSess,:),...
            'MarkerFaceAlpha', 0.5,...
            'MarkerEdgeAlpha', 0.5)
        axis square
        title(['Significant orientations: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
        xlabel([num2str(expt(1).experiment)])
        ylabel([num2str(expt(2).experiment)])
        refline(1)
        hold on 
    end
end

print(fullfile(fn_pop, ['sigOrisScatter_population.pdf']),'-dpdf','-bestfit')

%% get ODI

for iSess = 1:nSess
    for id = 1:nd
        [contra_resp{iSess,id} max_ind_contra{iSess,id}] = max(data_dfof_Ori{iSess,id}(:,:,find(Eyes),1),[],2);
        max_Ori_contra{iSess,id} = Oris(max_ind_contra{iSess,id});
        [ipsi_resp{iSess,id} max_ind_ipsi{iSess,id}] = max(data_dfof_Ori{iSess,id}(:,:,find(~Eyes),1),[],2); 
        max_Ori_ipsi{iSess,id} = Oris(max_ind_ipsi{iSess,id});

        real_contra_resp{iSess,id} = contra_resp{iSess,id};
        real_contra_ind{iSess,id} = find(real_contra_resp{iSess,id} < 0);
        real_contra_resp{iSess,id}(real_contra_ind{iSess,id}) = 0;

        real_ipsi_resp{iSess,id} = ipsi_resp{iSess,id};
        real_ipsi_ind{iSess,id} = find(real_ipsi_resp{iSess,id} < 0);
        real_ipsi_resp{iSess,id}(real_ipsi_ind{iSess,id}) = 0;

        ODI{iSess,id} = (real_contra_resp{iSess,id}-real_ipsi_resp{iSess,id})./(real_contra_resp{iSess,id}+real_ipsi_resp{iSess,id});
    end
end

%scatter of max response to contra and ipsi
figure
for id = 1:nd
    for iSess = 1:nSess
    subplot(2,2,id)
    contra_resp_any{iSess,id} = max(data_dfof_Ori{iSess,id}(resp_ind{iSess,id},:,find(Eyes),1),[],2);
    ipsi_resp_any{iSess,id} = max(data_dfof_Ori{iSess,id}(resp_ind{iSess,id},:,find(~Eyes),1),[],2);
%     C = linspace(1,10,length(contra_resp_any{iSess,id}(:,1)));
%     scatter(contra_resp_any{iSess,id},ipsi_resp_any{iSess,id},[],C)
    swarmchart(contra_resp_any{iSess,id},ipsi_resp_any{iSess,id}, sz_dots,...
            'MarkerFaceColor', color_choices(iSess,:),...
            'MarkerEdgeColor', color_choices(iSess,:),...
            'MarkerFaceAlpha', 0.5,...
            'MarkerEdgeAlpha', 0.5)
    hold on
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    axis square
    xlabel('Max df/f (Contra)')
    ylabel('Max df/f (Ipsi)')
    xlim([0.001 1])
    ylim([0.001 1])
    refline(1)
    title(['Responsive cells: ' num2str(expt(id).experiment)])
    end
end

% All ODI

for iSess = 1:nSess
    subplot(2,2,id+1)
    swarmchart(ODI{iSess,1}, ODI{iSess,2}, sz_dots,...
            'MarkerFaceColor', color_choices(iSess,:),...
            'MarkerEdgeColor', color_choices(iSess,:),...
            'MarkerFaceAlpha', 0.5,...
            'MarkerEdgeAlpha', 0.5)
    hold on
    axis square
    xlim([-1 1])
    ylim([-1 1])
    refline(1)
    xlabel([num2str(expt(1).experiment)])
    ylabel([num2str(expt(2).experiment)])
    title('ODI')
end

% Mean ODI
mean_ODI = zeros(nSess,nd);

for iSess = 1:nSess 
    for id = 1:nd
        mean_ODI(iSess,id) = nanmean(ODI{iSess,id});
    end
end

for iSess = 1:nSess
    subplot(2,2,id+2)
    swarmchart(mean_ODI(iSess,1), mean_ODI(iSess,2), sz_dots,...
            'MarkerFaceColor', color_choices(iSess,:),...
            'MarkerEdgeColor', color_choices(iSess,:),...
            'MarkerFaceAlpha', 0.5,...
            'MarkerEdgeAlpha', 0.5)
    hold on
    axis square
    xlim([-1 1])
    ylim([-1 1])
    refline(1)
    xlabel([num2str(expt(1).experiment)])
    ylabel([num2str(expt(2).experiment)])
    title('Mean ODI')
end

print(fullfile(fn_pop, ['OD_population.pdf']),'-dpdf','-bestfit')

save(fullfile(fn_pop, ['ODI_population.mat']), 'ODI', 'real_contra_resp', 'real_ipsi_resp', 'max_Ori_contra', 'max_Ori_ipsi')

%% Population averages for ODI

% % Do this for recovery, but histogram
% swarmchart(day_id(1:3),mean_ODI)
% axis square
% xlabel('Timepoint')
% set(gca,'XTick',3:5,'XTickLabel',{'Pre','Post','Rec'})
% ylabel('Mean ODI')
% ylim([-1 1])
% title('Mean ODI across days')

%% Population average for contra/ipsi inputs

% Mean contra/ipsi inputs
mean_ipsi = zeros(nSess,nd);
mean_contra = zeros(nSess,nd);

for iSess = 1:nSess
    for id = 1:nd
         mean_ipsi(iSess,id) = mean(real_ipsi_resp{iSess,id});
         mean_contra(iSess,id) = mean(real_contra_resp{iSess,id});
         mean_resp(iSess,id,:) = [mean_ipsi(iSess,id); mean_contra(iSess,id)]
    end
end

figure
for iSess = 1:nSess
    scatter([1 2 3 4],reshape(mean_resp(iSess,:,:),[1,4]), sz_dots,...
            'MarkerFaceColor', color_choices(iSess,:),...
            'MarkerEdgeColor', color_choices(iSess,:),...
            'MarkerFaceAlpha', 0.5,...
            'MarkerEdgeAlpha', 0.5)
    hold on
    axis square
    xlabel('Timepoint')
    set(gca,'XTick',1:4,'XTickLabel',{'Pre: Ipsi', 'Post: Ipsi','Pre: Contra','Post: Contra'})
    ylabel('Mean response (df/f)')
    title(['Ipsi v. contra response'])
    axis([0.5 4.5 0 0.35])
end

print(fullfile(fn_pop, ['EyeRespHist_population.pdf']),'-dpdf','-bestfit')

save(fullfile(fn_pop, ['EyeResp_population.mat']), 'mean_ipsi', 'mean_contra')

%% Need to finish editing this code
% %% von mises
% 
%     b_ori = cell(1,nd);
%     k1_ori = cell(1,nd);
%     R1_ori = cell(1,nd);
%     u1_ori = cell(1,nd);
%     R_square_ori = cell(1,nd);
%     sse_ori = cell(1,nd);
%     stim_OSI = cell(1,nd);
%     theta_hires = deg2rad(0:180);
%     y_fit = cell(1,nd);
% 
% for id = 1:nd    
%     for iEye = 1:nEye
%         for iCell = 1:nCells(iSess,id)
%             data = [data_dfof_Ori{id}(iCell,:,iEye) data_dfof_Ori{id}(iCell,1,iEye)];
%             theta = [deg2rad(Oris) pi];
%             [b_ori{id}(iEye,iCell),k1_ori{id}(iEye,iCell),R1_ori{id}(iEye,iCell),u1_ori{id}(iEye,iCell),sse_ori{id}(iEye,iCell),R_square_ori{id}(iEye,iCell)] ...
%                 = miaovonmisesfit_ori(theta,data);
%             [max_val{id} max_ind{id}] = max(data_dfof_Ori{id}(iCell,:,iEye),[],2);
%             null_ind{id} = max_ind{id}+(nOris./2);
%             null_ind{id}(find(null_ind{id}>nOris)) = null_ind{id}(find(null_ind{id}>nOris))-nOris;
%             min_val{id} = data_dfof_Ori{id}(iCell,null_ind{id},iEye);
%             if min_val{id}<0
%                 min_val{id} = 0;
%             end
%             stim_OSI{id}(1,iCell) = (max_val{id}-min_val{id})./(max_val{id}+min_val{id});
%             y_fit{id}(:,iEye,iCell) = b_ori{id}(iEye,iCell) + R1_ori{id}(iEye,iCell) .* exp(k1_ori{id}(iEye,iCell).*(cos(2.*(theta_hires-u1_ori{id}(iEye,iCell)))-1));
%         end
%     end
%     
%     [yfit_max{id}, yfit_max_ind{id}] = max(y_fit{id}(:,:,:),[],1);
%     prefOri_yfit{id} = squeeze(theta_hires(yfit_max_ind{id}));
% end   
% save(fullfile(fn_pop, ['oriResp_population.mat']), 'data_dfof_Ori','base_win','resp_win','h_Ori', 'b_ori', 'k1_ori', 'R1_ori', 'u1_ori', 'R_square_ori', 'sse_ori','stim_OSI', 'yfit_max', 'yfit_max_ind', 'prefOri_yfit') 
% 
% %% Population average for tuning.
% 
% % Plot tuning widths in cumulative distribution plot
% figure
% for iEye = 1:nEye
%     subplot(2,2,iEye)
%     for id = 1:nd
%         cdfplot(k1_ori{id}(iEye,:));
%         hold on
%     end
%     axis square
%     xlabel('Tuning width')
%     ylabel('Fraction of cells')
%     title(['Tuning width: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
% end
% legend('Pre','Post','Rec','Location','Southeast');
% 
% % Compare tuning widths of individual cells -- plot in scatterplots
% % (baseline v. MD)
% for iEye = 1:nEye
%     subplot(2,2,iEye+2)
%     scatter(k1_ori{1}(iEye,:),k1_ori{2}(iEye,:));
%     axis square
%     title(['Tuning width: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
%     xlabel([num2str(expt(1).experiment)])
%     ylabel([num2str(expt(2).experiment)])
%     refline(1)
%     hold on   
% end
% 
% print(fullfile(fn_pop, ['TuningPlots_population.pdf']),'-dpdf','-bestfit')
% 
% %% Compare dfof between pre- and post-MD
% 
% % EDIT: FINISH THIS SECTION
% figure
% for iEye = nEye
%     
% end
% 
% temp_data = data_dfof{1};
% temp_max = max(reshape(temp_data,[],size(temp_data,3)),[],1);
% 
% % EDIT: PUT TEMP_MAX IN HERE INSTEAD????
% [contra_resp{id} max_ind_contra{id}] = max(data_dfof_Ori{id}(:,:,find(Eyes),1),[],2);
% 
% subplot(2,1,iEye)
% scatter((data_dfof_max{1}(:,iEye)),(data_dfof_max{2}(:,iEye)))
% axis square
% title(['Max dfof: ' eye_str{find(contra(1:2)==Eyes(iEye))}]);
% xlabel([num2str(expt(1).experiment)])
% ylabel([num2str(expt(2).experiment)])
% refline(1)
>>>>>>> 56f5ff93 (Merge branch 'master' of github.com:Glickfeld-And-Hull-Laboratories/ImagingCode-Glickfeld-Hull)
% hold on   