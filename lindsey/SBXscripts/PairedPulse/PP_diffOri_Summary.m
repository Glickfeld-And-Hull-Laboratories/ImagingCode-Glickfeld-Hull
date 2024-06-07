
close all
clear all
clc

eval('PP_diffOri_2024_ExptList')
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey\Analysis\2P';
LG_base_out = fullfile(LG_base, 'Adaptation', 'ppDiffOri_2024');
if ~exist(LG_base_out)
    mkdir(LG_base_out)
end

%% collect datasets
nexp = size(expt,2);
mouse_mat = [];
max_dir_all = [];
pref_ori_all = [];
theta_90_all = [];
R_square_all = [];
max_dir_n_all = [];
fit_all = [];
OSI_all = [];
OSI_k_all = [];
HWHM_all = [];
delta_resp_all = [];
roc_resp_all = [];
pref_ori_sub_all = [];
fit_sub_all = [];
OSI_sub_all = [];
OSI_k_sub_all = [];
HWHM_sub_all = [];
delta_resp_sub_all = [];
roc_resp_sub_all = [];
nCells = [];
ppResp_all = cell(1,nexp);
ntot = zeros(1,nexp);
ori_sig_all = [];
ori_all = [];
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    mouse_mat = [mouse_mat mouse];
    date = expt(iexp).date;
    ImgFolder = expt(iexp).ImgFolder;
    nrun = size(ImgFolder,1);
    run_str = catRunName(ImgFolder, nrun);
    fprintf([mouse ' ' date '\n'])
    load(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_deltaResp.mat']))
    load(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respSig.mat']))
    max_dir_all = [max_dir_all; max_dir];
    pref_ori_all = [pref_ori_all; pref_ori];
    R_square_all = [R_square_all; R_square(:,:,1)];
    max_dir_n_all = [max_dir_n_all; max_dir_n'];
    ori_sig_all = [ori_sig_all; h2'];
    OSI_all = [OSI_all; OSI];
    OSI_k_all = [OSI_k_all; OSI_k];
    HWHM_all = [HWHM_all; HWHM];
    delta_resp_all = cat(1,delta_resp_all, delta_resp(:,:,:,1));
    fit_all = cat(2, fit_all, y_fit);
    pref_ori_sub_all = [pref_ori_sub_all; pref_ori_sub];
    OSI_sub_all = [OSI_sub_all; OSI_sub];
    OSI_k_sub_all = [OSI_k_sub_all; OSI_k_sub];
    HWHM_sub_all = [HWHM_sub_all; HWHM_sub];
    delta_resp_sub_all = cat(1,delta_resp_sub_all, delta_resp_sub(:,:,:,1));
    fit_sub_all = cat(2, fit_sub_all, y_fit_sub(:,:,:,1));
    theta_90_all = [theta_90_all; theta_90'];
    load(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_roc180v23.mat']))
    roc_resp_all = cat(1,roc_resp_all, roc_resp);
    roc_resp_sub_all = cat(1,roc_resp_sub_all, roc_resp_sub);
    ppResp_all{iexp} = ppResp;
    if ~exist(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_newFits.mat']))
        ori_fit = y_fit(:,:,3,1);
        save(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_newFits.mat']),'ori_fit', 'ppResp', 'R_square', 'theta_90')
    end
    ori_fit = y_fit(:,:,3,1);
    save(fullfile(LG_base_out, [date '_' mouse '_' run_str '_newFits.mat']),'ori_fit', 'ppResp', 'R_square', 'theta_90') 
end
load(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
load(fullfile(LG_base, [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
if length(ind_con)>0
    off_all = [offs; unique(cell2mat(input.tItiWaitFrames))];
    noff_all = length(off_all);
else
    off_all = offs;
    noff_all = noff;
end
control_ind = size(theta_90_all,2);
good_ind_theta = find(theta_90_all(:,control_ind)<= 22.5);
good_ind_dir = find(max_dir_n_all>= 500);
frameRateHz = double(frameRateHz);
mouse_list= mouse_mat;
save(fullfile(LG_base_out,'ppDiffOri_summary.mat'),'mouse_mat', 'max_dir_all', 'pref_ori_all', 'theta_90_all', 'OSI_all', 'delta_resp_all', 'fit_all', 'good_ind_theta');

%% summary of changes in tuning
% figure;
% col_mat = strvcat('b', 'r', 'y');
% start = 1;
% rep = 1;
% for iCell = 1:length(good_ind_theta)
%     iC = good_ind_theta(iCell);
%     if start>16
%         suptitle([mouse_list '- Cells ' num2str(iCell-16) ' to ' num2str(iCell-1) '- all fits'])
%         print(fullfile(LG_base_out, ['goodFits_summary' num2str(rep) '.pdf']),'-dpdf','-fillpage')
%         rep = 1+rep;
%         figure;
%         start = 1;
%     end
%     subplot(4, 4, start)
%     for ioff = 1:noff_all
%         plot(squeeze(fit_all(:,iC,ioff)),col_mat(ioff))
%         hold on
%         scatter(deltas, delta_resp_all(iC,:,ioff),['o' col_mat(ioff)])
%     end
%     start = start+1;
% end
% suptitle([mouse_list '- Cells ' num2str(iCell-start+2) ' to ' num2str(iCell) '- all fits'])
% print(fullfile(LG_base_out, ['goodFits_summary' num2str(rep) '.pdf']),'-dpdf','-fillpage')   

%% tuning figures
delta_resp_norm_all = bsxfun(@rdivide, delta_resp_all, max(delta_resp_all(:,:,3),[],2));
figure;
for idel = 1:nDelta
    ind_cells = intersect(good_ind_theta,find(max_dir_all == idel));
    for ioff = 1:noff_all
        subplot(noff_all,nDelta,((ioff-1)*nDelta)+idel)
        errorbar(deltas, mean(delta_resp_norm_all(ind_cells,:,ioff),1), std(delta_resp_norm_all(ind_cells,:,ioff),[],1)./sqrt(length(ind_cells)), 'ok')
        if ioff == 1
            title([num2str(deltas(idel)) ' deg pref- ' num2str(length(ind_cells)) ' cells'])
        end
        ylim([-0.1 1.2])
    end
end
suptitle([mouse_list '- avg all cells- normalized'])
print(fullfile(LG_base_out, 'tuning_byPref_byInt_sep_summary.pdf'),'-dpdf','-fillpage')

figure;
c = [.6 .6 .6; .4 .4 .4; 0 0 0];
[n n2] = subplotn(nDelta);
for idel = 1:nDelta
    subplot(n,n2,idel)
    ind_cells = intersect(good_ind_theta,find(max_dir_all == idel));
    for ioff = 1:noff_all
        errorbar(deltas, mean(delta_resp_norm_all(ind_cells,:,ioff),1), std(delta_resp_norm_all(ind_cells,:,ioff),[],1)./sqrt(length(ind_cells)), 'o', 'Color', c(ioff,:));
        hold on
        if ioff == 1
            title([num2str(deltas(idel)) ' deg pref- ' num2str(length(ind_cells)) ' cells'])
        end
        ylim([-0.1 1.2])
    end
    xlabel('Stimulus Orientation (deg)')
end
suptitle([mouse_list '- avg all cells- normalized'])
print(fullfile(LG_base_out, 'tuning_byPref_byInt_summary.pdf'),'-dpdf','-fillpage')

figure;
c = [.6 .6 .6; .4 .4 .4; 0 0 0];
[n n2] = subplotn(nDelta);
for idelta = 1:nDelta
    for ioff = 1:noff_all
        for idel = 1:nDelta
            ind_cells = intersect(good_ind_theta, find(max_dir_all == idel));
            subplot(n,n2,idelta)
            errorbar(deltas(idel), mean(delta_resp_norm_all(ind_cells,idelta,ioff),1),std(delta_resp_norm_all(ind_cells,idelta,ioff),[],1)./sqrt(length(ind_cells)),'o', 'Color', c(ioff,:));
            hold on
        end
    end
    ylim([-0.1 1.2])
    xlabel('Cell preference group (deg)')
    title([num2str(deltas(idelta)) ' deg change'])
end
suptitle([mouse_list '- avg all cells- normalized'])
print(fullfile(LG_base_out, 'popTuning_byInt_summary.pdf'),'-dpdf','-fillpage')


%% Pref and OSI change
deltas = 22.5:22.5:180;
delta_diff = 180-deltas;
delta_diff(find(delta_diff>90)) = 180- delta_diff(find(delta_diff>90));
diffs = unique(delta_diff);
ndiff = length(diffs);
pref_ori_all_diff = pref_ori_all;
pref_ori_all_diff(find(pref_ori_all_diff>90)) = 180-pref_ori_all_diff(find(pref_ori_all_diff>90));


cell_group_n = zeros(1,ndiff);
group_list = zeros(size(max_dir_all,1),1);
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    cell_group_n(1,idiff) = length(cell_ind);
    group_list(cell_ind,1) = idiff;
end

col_mat = strvcat('b', 'r', 'y');
figure;
subplot(2,2,1)
osi_diff = zeros(noff, ndiff,2);
osi_avg = zeros(ndiff,2);
p_osi_diff = zeros(noff,ndiff);
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    osi_avg(idiff,1) = mean(OSI_all(cell_ind,3),1);
    osi_avg(idiff,2) = std(OSI_all(cell_ind,3),[],1)./sqrt(length(cell_ind));
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_all(cell_ind,ioff)-OSI_all(cell_ind,3),1), std(OSI_all(cell_ind,ioff)-OSI_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
        osi_diff(ioff,idiff,1) = mean(OSI_all(cell_ind,ioff)-OSI_all(cell_ind,3),1);
        osi_diff(ioff,idiff,2) = std(OSI_all(cell_ind,ioff)-OSI_all(cell_ind,3),[],1)./sqrt(length(cell_ind));
        [h, p_osi_diff(ioff,idiff)] = ttest(OSI_all(cell_ind,ioff)-OSI_all(cell_ind,3));
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in OSI'])
title('OSI')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
subplot(2,2,2)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_all(cell_ind,ioff)./OSI_all(cell_ind,3),1), std(OSI_all(cell_ind,ioff)./OSI_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Ratio of OSI'])
title('OSI')
ylim([0.5 1.5])
xlim([-10 100])
hline(1)
subplot(2,2,3)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_k_all(cell_ind,ioff)-OSI_k_all(cell_ind,3),1), std(OSI_k_all(cell_ind,ioff)-OSI_k_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in OSI-k'])
title('OSI-k')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
subplot(2,2,4)
HWHM_deg_all = rad2deg(HWHM_all);
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(HWHM_deg_all(cell_ind,ioff)-HWHM_deg_all(cell_ind,3),1), std(HWHM_deg_all(cell_ind,ioff)-HWHM_deg_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in HWHM'])
title('HWHM')
ylim([-30 30])
xlim([-10 100])
hline(0)
suptitle([mouse_list '- OSI by Interval'])
print(fullfile(LG_base_out, 'OSI_byInt_summary.pdf'),'-dpdf','-fillpage')

figure;
subplot(2,2,1)
g = jet(180);
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    plot(pref_ori_all(iC, 2),pref_ori_all(iC, 1), 'o', 'Color', g(ceil(pref_ori_all(iC,3))+1,:))
    hold on
end
refline(1,0)
axis square
xlabel(['Pref Ori- ' num2str(chop(off_all(2).*(1000/frameRateHz),3)) ' ms ISI'])
ylabel(['Pref Ori- ' num2str(chop(off_all(1).*(1000/frameRateHz),3)) ' ms ISI'])
subplot(2,2,2)
g = jet(91);
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    plot(pref_ori_all_diff(iC, 2),pref_ori_all_diff(iC, 1), 'o', 'Color', g(ceil(pref_ori_all_diff(iC,3)+1),:))
    hold on
end
refline(1,0)
axis square
xlabel(['Diff Pref from Adapt Ori- ' num2str(chop(off_all(2).*(1000/frameRateHz),3)) ' ms ISI'])
ylabel(['Diff Pref from Adapt Ori- ' num2str(chop(off_all(1).*(1000/frameRateHz),3)) ' ms ISI'])
pref_ori_diff = zeros(noff, ndiff,2);
p_pref_ori_diff = zeros(noff,ndiff);
subplot(2,2,3)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),1), std(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
        pref_ori_diff(ioff,idiff,1) = mean(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),1);
        pref_ori_diff(ioff,idiff,2) = std(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),[],1)./sqrt(length(cell_ind));
        [h, p_pref_ori_diff(ioff,idiff)] = ttest(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3));
    end
end
xlabel(['Diff of Max from Adaptor (deg)'])
ylabel(['Change in Pref'])
ylim([-40 40])
hline(0)
xlim([-10 100])
suptitle([mouse_list '- Change in Preference by Interval'])
print(fullfile(LG_base_out, 'prefDiff_byInt_summary.pdf'),'-dpdf','-fillpage')

col_mat = strvcat('b','r','y');
figure;
delta_resp_all_avg = zeros(noff_all, ndiff,2);
delta_resp_norm_all_avg = zeros(noff_all, ndiff,2);
p_delta_resp_norm = zeros(noff, ndiff);
for i = 1:noff_all
    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        subplot(2,2,1)
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        delta_resp_all_avg(i,idiff,1) = squeeze(mean(mean(delta_resp_all(good_ind_theta,del,i),2),1));
        delta_resp_all_avg(i,idiff,2) = squeeze(std(mean(delta_resp_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta));
        hold on
        subplot(2,2,2)
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_norm_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_norm_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        hold on
        delta_resp_norm_all_avg(i,idiff,1) = squeeze(mean(mean(delta_resp_norm_all(good_ind_theta,del,i),2),1));
        delta_resp_norm_all_avg(i,idiff,2) = squeeze(std(mean(delta_resp_norm_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta));
        if i < 3
            [h p_delta_resp_norm(i,idiff)] = ttest(mean(delta_resp_norm_all(good_ind_theta,del,3),2), mean(delta_resp_norm_all(good_ind_theta,del,i),2));
        end
    end
end
subplot(2,2,1)
title('Absolute')
ylabel('Mean dF/F')
xlabel('Diff of Stim from Adaptor (deg)')
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([0 .3])
subplot(2,2,2)
title('Normalized')
ylabel('Mean dF/F')
xlabel('Diff of Stim from Adaptor (deg)')
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([0 1])
subplot(2,2,3)
delta_resp_group_avg = zeros(noff, ndiff,2);
p_delta_resp_group = zeros(noff,ndiff);
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    diff_resp_group= [];
    for i = 1:length(del)
        cell_ind = intersect(good_ind_theta, find(max_dir_all == del(i)));
        diff_resp_group = [diff_resp_group; squeeze(delta_resp_norm_all(cell_ind,del(i),:))];
    end
    for i = 1:noff
        errorbar(diffs(idiff), mean(diff_resp_group(:,i),1), std(diff_resp_group(:,i),[],1)./sqrt(size(diff_resp_group,1)),['o' col_mat(i,:)])
        hold on
        delta_resp_group_avg(i,idiff,1) = mean(diff_resp_group(:,i),1);
        delta_resp_group_avg(i,idiff,2) = std(diff_resp_group(:,i),[],1)./sqrt(size(diff_resp_group,1));
        [h p_delta_resp_group(i,idiff)] = ttest(diff_resp_group(:,3),diff_resp_group(:,i));
    end
end
title('Normalized')
ylabel('dF/F at max ori')
xlabel('Diff of Max from Adaptor (deg)')
ylim([0 2])
xlim([-10 100])
hline(1)
suptitle([mouse_list '- Mean Resp by Interval'])
print(fullfile(LG_base_out, 'resp_byInt_summary.pdf'),'-dpdf','-fillpage')

figure;
subplot(2,2,1)
nCells = length(good_ind_theta);
delta_resp_sum = zeros(noff*nCells,ndiff);
start = 1;
off_i = [];
diff_i = [];
for i = 1:noff
    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_all(good_ind_theta,del,3),2)-mean(delta_resp_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_all(good_ind_theta,del,3),2)-mean(delta_resp_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        delta_resp_sum(start:start+nCells-1,idiff) = mean(delta_resp_all(good_ind_theta,del,3),2)-mean(delta_resp_all(good_ind_theta,del,i),2);
        hold on
    end
    start = start+nCells;
end
ylabel('Diff from control resp')
xlabel('Diff of Stim from Adaptor (deg)')
title('Abs')
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([0 0.2])
subplot(2,2,2)
delta_resp_norm_sum = zeros(noff*nCells,ndiff);
diff_resp_norm = zeros(size(delta_resp_norm_all,1),ndiff,noff_all);
start = 1;
off_ii = [];
diff_ii = [];
for i = 1:noff
    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_norm_all(good_ind_theta,del,3),2)-mean(delta_resp_norm_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_norm_all(good_ind_theta,del,3),2)-mean(delta_resp_norm_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        diff_resp_norm(:,idiff,i) = mean(delta_resp_norm_all(:,del,i),2);
        if i == 1
            diff_resp_norm(:,idiff,3) = mean(delta_resp_norm_all(:,del,3),2);
        end
        hold on
        delta_resp_norm_sum(start:start+nCells-1,idiff) = mean(delta_resp_norm_all(good_ind_theta,del,3),2)-mean(delta_resp_norm_all(good_ind_theta,del,i),2);
    end
    start = start+nCells;
end
ylabel('Diff from control resp')
xlabel('Diff of Stim from Adaptor (deg)')
title('Norm')
hline(0)
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([-0.2 0.4])

[n, bin] = histc(pref_ori_all_diff(:,3),[0 15 75 91]);
subplot(2,2,3)
delta_resp_group = nan(sum(noff.*n,1),1);
start = 1;
off_id = [];
diff_id = [];
for idiff = 1:length(n)-1
    ind = intersect(good_ind_theta, find(bin == idiff));
    diff_resp_group= zeros(length(ind),noff);
    for i = 1:length(ind)
        [peak_val peak_loc] = max(fit_all(:,ind(i),3),[],1);
        fit_norm = squeeze(fit_all(:,ind(i),:))./peak_val;
        diff_resp_group(i,:) = fit_norm(peak_loc,1:2);
    end
    for i = 1:noff
        errorbarxy(mean(pref_ori_all_diff(ind,3),1), mean(diff_resp_group(:,i),1), std(pref_ori_all_diff(ind,3),[],1)./sqrt(size(diff_resp_group,1)),std(diff_resp_group(:,i),[],1)./sqrt(size(diff_resp_group,1)),{['o' col_mat(i,:)], col_mat(i,:),col_mat(i,:)})
        hold on
        delta_resp_group(start:start+length(ind)-1,1) = diff_resp_group(:,i);
        off_id = [off_id ones(1,length(ind))*offs(i)];
        diff_id = [diff_id ones(1,length(ind))*diffs(idiff)];
        start = start+length(ind);
    end
end
title('Norm')
ylabel('dF/F at max ori')
xlabel('Diff of Pref from Adaptor (deg)')
ylim([0 2])
xlim([-10 100])
hline(1)
suptitle([mouse_list '- Diff Resp by Stimulus'])
print(fullfile(LG_base_out, 'diffResp_byStim_byFit_newBin2_summary.pdf'),'-dpdf','-fillpage')

[p_sum, table_sum, stats_sum] = anova2(delta_resp_sum, nCells);
[p_normsum, table_normsum, stats_normsum] = anova2(delta_resp_norm_sum, nCells);
[p_normgroup, table_normgroup, stats_normgroup] = anovan(delta_resp_group, {off_id, diff_id});

figure;
subplot(2,2,1)
pref_ori_group = nan(sum(noff.*n,1),1);
start = 1;
for idiff = 1:length(n)-1
    ind = intersect(good_ind_theta, find(bin == idiff));
    for ioff = 1:noff
        errorbarxy(mean(pref_ori_all_diff(ind,3),1), mean(pref_ori_all_diff(ind,ioff)-pref_ori_all_diff(ind,3),1), std(pref_ori_all_diff(ind,3),[],1)./sqrt(length(ind)), std(pref_ori_all_diff(ind,ioff)-pref_ori_all_diff(ind,3),[],1)./sqrt(length(ind)), {['o' col_mat(ioff,:)], col_mat(ioff,:),col_mat(ioff,:)})
        hold on
        pref_ori_group(start:start+length(ind)-1,1) = pref_ori_all_diff(ind,ioff)-pref_ori_all_diff(ind,3);
        start= start +length(ind);
    end
end
xlabel(['Diff of Pref from Adaptor (deg)'])
ylabel(['Difference in Pref (deg)'])
title('Pref')
ylim([-20 50])
hline(0)
xlim([-10 100])
subplot(2,2,2)
OSI_group = nan(sum(noff.*n,1),1);
start = 1;
for idiff = 1:length(n)-1
    ind = intersect(good_ind_theta, find(bin == idiff));
    for ioff = 1:noff
        errorbarxy(mean(pref_ori_all_diff(ind,3),1), mean(OSI_all(ind,ioff)-OSI_all(ind,3),1),std(pref_ori_all_diff(ind,3),[],1)./sqrt(length(ind)), std(OSI_all(ind,ioff)-OSI_all(ind,3),[],1)./sqrt(length(ind)), {['o' col_mat(ioff,:)], col_mat(ioff,:),col_mat(ioff,:)})
        hold on
        OSI_group(start:start+length(ind)-1,1) = OSI_all(ind,ioff)-OSI_all(ind,3);
        start= start +length(ind);
    end
end
xlabel(['Diff of Pref from adaptor (deg)'])
ylabel(['Difference in OSI'])
title('OSI')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
subplot(2,2,3)
k = nan(sum(noff.*n,1),1);
start = 1;
for idiff = 1:length(n)-1
    ind = intersect(good_ind_theta, find(bin == idiff));
    for ioff = 1:noff
        errorbarxy(mean(pref_ori_all_diff(ind,3),1), mean(OSI_k_all(ind,ioff)-OSI_k_all(ind,3),1),std(pref_ori_all_diff(ind,3),[],1)./sqrt(length(ind)), std(OSI_k_all(ind,ioff)-OSI_k_all(ind,3),[],1)./sqrt(length(ind)), {['o' col_mat(ioff,:)], col_mat(ioff,:),col_mat(ioff,:)})
        hold on
        OSIk_group(start:start+length(ind)-1,1) = OSI_k_all(ind,ioff)-OSI_k_all(ind,3);
        start= start +length(ind);
    end
end
xlabel(['Diff of Pref from adaptor (deg)'])
ylabel(['Difference in OSI-k'])
title('OSI-k')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
suptitle([mouse_list '- Change in Preference by Interval'])
print(fullfile(LG_base_out, 'diff_byInt_byPref_newBin_summary.pdf'),'-dpdf','-fillpage')
[p_prefdiff, table_prefdiff, stats_prefdiff] = anovan(pref_ori_group, {off_id, diff_id});
[p_osidiff, table_osidiff, stats_osidiff] = anovan(OSI_group, {off_id, diff_id});

save(fullfile(LG_base_out, 'tuningStatsSummary.mat'),'table_osidiff','table_prefdiff','table_normsum','table_normgroup')
%% max log likelihood

nboot = 1000;
OSI_cells = good_ind_theta;
cellN = length(OSI_cells);
sz_fit = size(fit_all,1);

loglike_all_fun = nan(sz_fit,nDelta,noff_all,nboot+1);
loglike_all_fun_nosub = nan(sz_fit,nDelta,noff_all,nboot+1);
loglike_all_fun_subfact = nan(sz_fit,nDelta,noff_all,nboot+1);
maxloglike_all = nan(nDelta,noff_all,nboot+1);
maxloglike_all_subfact = nan(nDelta,noff_all,nboot+1);
maxloglike_all_nosub = nan(nDelta,noff_all,nboot+1);

[max_fit, max_ind] = max(fit_all(:,:,3),[],1);
delta_resp_all_thresh = bsxfun(@rdivide,delta_resp_all,max_fit').*100;
delta_resp_all_thresh(find(delta_resp_all_thresh<1)) = 1;
delta_resp_all_thresh(find(delta_resp_all_thresh>170)) = 170;
fit_all_thresh =  bsxfun(@rdivide,fit_all(:,:,3),max_fit)*100;
fit_all_thresh(find(fit_all_thresh<1)) = 1;
fit_all_thresh(find(fit_all_thresh>170)) = 170;
% [n_fit bin_fit] = histc(pref_ori_all(:,3), [0:22.5:180]);
% delta_avg = zeros(length(n_fit)-1, nDelta,noff_all);
% fit_avg = zeros(sz_fit, length(n_fit)-1);
% delta_avg_thresh = zeros(length(n_fit)-1, nDelta,noff_all);
% fit_avg_thresh = zeros(sz_fit, length(n_fit)-1);
% loglike_all_neurons = nan(length(n_fit)-1,sz_fit,nDelta,noff_all);
% loglike_all_fact = nan(length(n_fit)-1,1,nDelta,noff_all);
% loglike_all_sum = nan(length(n_fit)-1,sz_fit,nDelta,noff_all);
% for i = 1:length(n_fit)-1
%     ind = find(bin_fit == i);
%     fit_avg(:,i) = squeeze(mean(fit_all(:,ind,3),2))';
%     delta_avg(i,:,:) = mean(delta_resp_all(ind,:,:),1);
%     max_fit = max(fit_avg(:,i),[],1);
%     fit_avg_thresh(:,i) = (fit_avg(:,i)./max_fit)*100;
%     fit_avg_thresh(find(fit_avg_thresh<1)) = 1;
%     fit_avg_thresh(find(fit_avg_thresh>170)) = 170;
%     delta_avg_thresh(i,:,:) = (delta_avg(i,:,:)./max_fit)*100;
%     delta_avg_thresh(find(delta_avg_thresh<1)) = 1;
%     delta_avg_thresh(find(delta_avg_thresh>170)) = 170;
%     for idelta = 1:nDelta
%         for ioff =1:noff_all
%             loglike_all_neurons(i,:,idelta,ioff) = squeeze(log10(fit_avg_thresh(:,i)').*delta_avg_thresh(i,idelta,ioff))';
%             loglike_all_sum(i,:,idelta,ioff) = squeeze(fit_avg_thresh(:,i)');
%             loglike_all_fact(i,1,idelta,ioff) = log10(gamma(delta_avg_thresh(i,idelta,ioff)+1));
%         end
%     end
% end
% loglike_all_fun =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons,1), bsxfun(@plus,(nansum(loglike_all_sum,1)),nansum(loglike_all_fact,1))));
% loglike_all_fun_nosub =  squeeze(nansum(loglike_all_neurons,1));
% [max_val maxloglike_all] = max(loglike_all_fun,[],1);
% [max_val maxloglike_all_nosub] = max(loglike_all_fun_nosub,[],1);
% maxloglike_all = squeeze(maxloglike_all);
% maxloglike_all_nosub = squeeze(maxloglike_all_nosub);
loglike_all_neurons = nan(cellN,sz_fit,nDelta,noff_all);
loglike_all_fact = nan(cellN,1,nDelta,noff_all);
loglike_all_sum = nan(cellN,sz_fit,nDelta,noff_all);
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    for idelta = 1:nDelta
        for ioff =1:noff_all
            loglike_all_neurons(iC,:,idelta,ioff) = squeeze(log10(fit_all_thresh(:,iC)').*delta_resp_all_thresh(iC,idelta,ioff))';
            loglike_all_sum(iC,:,idelta,ioff) = squeeze(fit_all_thresh(:,iC)');
            loglike_all_fact(iC,1,idelta,ioff) = log10(gamma(delta_resp_all_thresh(iC,idelta,ioff)+1));
        end
    end
end

%find best subtraction to account for inhomogeneity
loglike_all_fish = nan(sz_fit,nDelta,100);
maxloglike_fish = nan(nDelta,100);
for n = 1:100;
    loglike_all_fish(:,:,:,n) =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons(OSI_cells,:,:,:),1), bsxfun(@plus,(nansum(loglike_all_sum(OSI_cells,:,:,:),1)./n),nansum(loglike_all_fact(OSI_cells,:,:,:),1))));
    for idelta = 1:nDelta
        [max_val maxloglike_fish(idelta,n)] = max(loglike_all_fish(:,idelta,3,n),[],1);
    end
end
maxloglike_temp = maxloglike_fish;
ind = find(maxloglike_temp(8,:)<90);
maxloglike_temp(8,ind) = maxloglike_temp(8,ind)+180;
maxloglike_diff = bsxfun(@minus,maxloglike_temp, deltas');
maxloglike_sos = sum(maxloglike_diff.^2,1);
[min_val, min_ind] = min(maxloglike_sos,[],2);

for iboot = 1:nboot+1
    if iboot>1
        ind_cells_temp = OSI_cells(randsample(cellN, cellN, 1));
    else
        ind_cells_temp = OSI_cells;
    end
    loglike_all_fun(:,:,:,iboot) =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons(ind_cells_temp,:,:,:),1), bsxfun(@plus,(nansum(loglike_all_sum(ind_cells_temp,:,:,:),1)./min_ind),nansum(loglike_all_fact(ind_cells_temp,:,:,:),1))));
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            [max_val maxloglike_all(idelta,ioff,iboot)] = max(loglike_all_fun(:,idelta,ioff,iboot),[],1);
        end
    end
end

% figure;
% for idelta = 1:nDelta
%     subplot(3,3,idelta)
%     for ioff = 1:noff_all
%         plot([0:180], loglike_all_fun(:,idelta,ioff));
%         hold on
%     end
%     title([num2str(deltas(idelta)) 'deg'])
%     xlabel('Orientation')
%     ylabel('Log likelihood')
%     xlim([-10 190])
%     %ylim([-3000 0])
% end

figure;
for idelta = 1:nDelta
    subplot(3,3,idelta)
    for ioff = 1:noff_all
        shadedErrorBar([0:180], loglike_all_fun(:,idelta,ioff,1), squeeze(std(loglike_all_fun(:,idelta,ioff,2:end),[],4)), {[col_mat(ioff) '-'],'markerfacecolor',col_mat(ioff)});
        hold on
    end
    title([num2str(deltas(idelta)) 'deg'])
    xlabel('Orientation')
    ylabel('Log likelihood')
    xlim([-10 190])
    %ylim([-3000 0])
end
suptitle([mouse_list '- Log Likelihood Function- all fits- 100X'])
print(fullfile(LG_base_out, 'loglike_byDelta_byInt_summary_allFits.pdf'),'-dpdf','-fillpage')



%likelihood is 0
figure;
for ioff = 1:noff_all
    temp_likely = squeeze(loglike_all_fun(:,:,ioff,1));
    errorbar([0 deltas], [temp_likely(end,end) temp_likely(end,:)],[squeeze(std(loglike_all_fun(end,end,ioff,2:end),[],4)) squeeze(std(loglike_all_fun(end,:,ioff,2:end),[],4))]);
    hold on
    xlim([-22 202])
    xlabel('Orientation')
    ylabel('Log likelihood of 0 deg')
end
suptitle([mouse_list '- Log Likelihood of 0 deg stimulus- all cells- 100X'])
print(fullfile(LG_base_out, 'loglike_is0_combo_allCells_100X.pdf'),'-dpdf','-fillpage')

figure;
for ioff = 1:noff_all
    subplot(2,1,1)
    plot([0 deltas], [maxloglike_all(end,ioff,1);maxloglike_all(:,ioff,1)],'-')
    hold on
    subplot(2,1,2)
    scatter(reshape(repmat([0 deltas], [nboot, 1])', [1, nboot*(nDelta+1)]), reshape(squeeze([maxloglike_all(end,ioff,2:nboot+1);maxloglike_all(:,ioff,2:nboot+1)]), [1, nboot*(nDelta+1)]),'o')
    hold on
end
for i = 1:2
    subplot(2,1,i)
    xlabel('Orientation')
    ylabel('Max log likelihood')
    axis square
end
subplot(2,1,1)
title('Orig')
subplot(2,1,2)
title('Bootstrap')
suptitle([mouse_list '- Max Log Likelihood - scaled sub- 100X'])
print(fullfile(LG_base_out, 'maxloglike_scaledsub_100X.pdf'),'-dpdf','-fillpage')

maxloglike_change = maxloglike_all;
maxloglike_change(find(maxloglike_all>90)) = 179-maxloglike_change(find(maxloglike_all>90));
maxloglike_change_avg = zeros(ndiff,noff_all,2);
maxloglike_change_all = zeros(ndiff*nboot, noff);
start = 1;
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    maxloglike_change_avg(idiff,:,1) = mean(maxloglike_change(del,:,1),1);
    maxloglike_change_avg(idiff,:,2) = std(mean(maxloglike_change(del,:,2:nboot+1),1),[],3);
    maxloglike_change_all(start:start+nboot-1,:) = squeeze(mean(maxloglike_change(del,1:2,2:nboot+1),1))';
    start = start+nboot;
end
figure;
for i = 1:noff_all
	errorbar(diffs, maxloglike_change_avg(:,i,1), maxloglike_change_avg(:,i,2),'-o');
    hold on
end
xlabel('Stim - Adapter')
ylabel('Max likely ori - Adapter')
xlim([-10 100])
ylim([-10 100])
suptitle([mouse_list '- Max Log Likelihood dist from Adapter - scaled sub- 100X'])
print(fullfile(LG_base_out, 'maxloglike_distAdapt_scaledsub_100X.pdf'),'-dpdf','-fillpage')
[p_max, table_max, stats_max] = anova2(maxloglike_change_all, nboot);


figure;
delta_sq_boot = zeros(nDelta,nDelta,noff_all);
delta_sq = zeros(nDelta,nDelta,noff_all);
for ioff = 1:noff_all
    for idelta = 1:nDelta
        for idel = 1:nDelta
            delta_sq_boot(idel,idelta,ioff) = length(find(maxloglike_all(idelta,ioff,2:end) == idel))./1000;
        end
        delta_sq(maxloglike_all(idelta,ioff,1),idelta,ioff) = 1;
    end
    subplot(3,2,(ioff*2)-1)
    imagesc(flipud(delta_sq(:,:,ioff)))
    axis square
    set(gca, 'XTick', 1:nDelta, 'XTickLabels',deltas)
    set(gca, 'YTick', 1:nDelta, 'YTickLabels',fliplr(deltas))
    xlabel('Actual Ori (deg)')
    ylabel('Max Likelihood Ori (deg)')
    title([num2str(chop(off_all(ioff).*(1000/frameRateHz),3)) ' ms ISI'])
    colormap(hot)
    subplot(3,2,(ioff*2))
    imagesc(flipud(delta_sq_boot(:,:,ioff)),[0 1])
    axis square
    set(gca, 'XTick', 1:nDelta, 'XTickLabels',deltas)
    set(gca, 'YTick', 1:nDelta, 'YTickLabels',fliplr(deltas))
    title([num2str(chop(off_all(ioff).*(1000/frameRateHz),3)) ' ms ISI- bootstrap'])
    xlabel('Actual Ori (deg)')
    ylabel('Max Likelihood Ori (deg)')
    colormap(hot)
    colorbar
end
suptitle([mouse_list '- Max Log Likelihood- All cells- 100X'])
print(fullfile(LG_base_out, 'maxloglike_byInt_summary_allCells.pdf'),'-dpdf','-fillpage')

save(fullfile(LG_base_out, 'maxLogLike_summary.mat'),'delta_sq_boot','delta_sq','maxloglike_all','loglike_all_neurons','loglike_all_sum', 'loglike_all_fact')

%% max log likelihood by trials
col_mat = strvcat('b','r','y');
OSI_cells = good_ind_theta;
cellN = length(OSI_cells);
sz_fit = size(fit_all,1);
n_ind = zeros(nexp, nDelta,noff_all);
for iexp = 1:nexp
    ppTemp = ppResp_all{iexp};
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            n_ind(iexp,idelta,ioff) = size(ppTemp{ioff, idelta},2);
        end
    end
end
min_trials = squeeze(min(n_ind,[],1));

[max_fit, max_ind] = max(fit_all(:,:,3),[],1);
fit_all_thresh =  bsxfun(@rdivide,fit_all(:,:,3),max_fit)*100;
fit_all_thresh(find(fit_all_thresh<1)) = 1;
fit_all_thresh(find(fit_all_thresh>170)) = 170;

ppResp_trials = cell(noff_all,nDelta);
start = 1;
for iexp = 1:nexp
    ppTemp = ppResp_all{iexp};
    nc = size(ppTemp{1,1},1);
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            ppT = ppTemp{ioff,idelta};
            ppT_norm = bsxfun(@rdivide, ppT(:,1:min_trials(idelta,ioff)), max_fit(1,start:nc+start-1)')*100;
            ppT_norm(find(ppT_norm<1)) = 1;
            ppT_norm(find(ppT_norm>170)) = 170;
            ppResp_trials{ioff,idelta}(start:nc+start-1,:) = ppT_norm;
        end
    end
    start = start+nc;
end

cellN = length(good_ind_theta);
nCells = size(fit_all_thresh,2);
nboot = 100;
sz_fit = size(fit_all,1);
loglike_all_neurons = cell(nDelta,noff_all);
loglike_all_fact = cell(nDelta,noff_all);
loglike_all_sum = cell(nDelta,noff_all);
loglike_all_fun = cell(nDelta,noff_all,nboot+1);
maxloglike_all = cell(nDelta,noff_all,nboot+1);
maxloglike_change = cell(nDelta,noff_all,nboot+1);

for idelta = 1:nDelta
    for ioff = 1:noff_all
        loglike_all_neurons{idelta,ioff} = zeros(sz_fit,nCells,min_trials(idelta,ioff));
        loglike_all_sum{idelta,ioff} = zeros(sz_fit,nCells,min_trials(idelta,ioff));
        loglike_all_fact{idelta,ioff} = zeros(1,nCells,min_trials(idelta,ioff));
        for i = 1:min_trials(idelta,ioff)
            for iCell = 1:cellN
                iC = good_ind_theta(iCell);
                loglike_all_neurons{idelta,ioff}(:,iC,i) = squeeze(bsxfun(@times,log(fit_all_thresh(:,iC)'),ppResp_trials{ioff,idelta}(iC,i)))';
                loglike_all_sum{idelta,ioff}(:,iC,i) = fit_all_thresh(:,iC);
                loglike_all_fact{idelta,ioff}(1,iC,i) = log(gamma(ppResp_trials{ioff,idelta}(iC,i)+1));
            end
        end
        for iboot = 1:nboot+1
            fprintf('.')
            if iboot>1
                ind_cells_temp = good_ind_theta(randsample(cellN, cellN, 1));
            else
                ind_cells_temp = good_ind_theta;
            end
   
            loglike_all_fun{idelta,ioff,iboot} =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons{idelta,ioff}(:,ind_cells_temp,:),2), bsxfun(@plus,(nansum(loglike_all_sum{idelta,ioff}(:,ind_cells_temp,:),2)),nansum(loglike_all_fact{idelta,ioff}(:,ind_cells_temp,:),2))));
            [max_val maxloglike_all{idelta,ioff,iboot}] = max(loglike_all_fun{idelta,ioff,iboot},[],1);
            maxloglike_change{idelta,ioff,iboot} = maxloglike_all{idelta,ioff,iboot}-1;
            maxloglike_change{idelta,ioff,iboot}(find(maxloglike_all{idelta,ioff,iboot}>90)) = 179-maxloglike_change{idelta,ioff,iboot}(find(maxloglike_all{idelta,ioff,iboot}>90));
        end
    end

    if iboot == 0
        figure;
        for idelta = 1:nDelta
            for ioff = 1:noff_all
                subplot(3,3,idelta)
                shadedErrorBar(0:180, mean(loglike_all_fun{idelta,ioff,iboot},2), std(loglike_all_fun{idelta,ioff,iboot},[],2)./sqrt(size(loglike_all_fun{idelta,ioff,iboot},2)),{[col_mat(ioff) '-'],'markerfacecolor',col_mat(ioff)});
                hold on
            end
            title([num2str(deltas(idelta)) 'deg'])
            xlabel('Orientation')
            ylabel('Log likelihood')
            xlim([-10 190])
            %ylim([-3000 0])
        end
        suptitle([mouse_list '- Log Likelihood Function- all fits- by trial- 100X'])
        print(fullfile(LG_base_out, 'loglike_byDelta_byInt_summary_allFits_byTrial.pdf'),'-dpdf','-fillpage')

        %likelihood is 0
        figure;
        likely_all = [];
        off_id = [];
        for ioff = 1:noff_all
            for idelta = 1:nDelta
                temp_likely = loglike_all_fun{idelta,ioff,iboot};
                errorbar(deltas(idelta), mean(temp_likely(end,:),2),std(temp_likely(end,:),[],2)./sqrt(size(temp_likely,2)),['o' col_mat(ioff,:)]);
                hold on
                if idelta == nDelta
                    errorbar(0, mean(temp_likely(end,:),2),std(temp_likely(end,:),[],2)./sqrt(size(temp_likely,2)),['o' col_mat(ioff,:)]);
                end
                if idelta == 1
                    likely_all = [likely_all temp_likely(end,:)];
                    off_id = [off_id ioff*ones(size(temp_likely(1,:)))];
                end
            end
            xlim([-22 202])
            xlabel('Orientation')
            ylabel('Log likelihood of 0 deg')
        end
        suptitle([mouse_list '- Log Likelihood of 0 deg stimulus- all cells- by trial- 100X'])
        print(fullfile(LG_base_out, 'loglike_is0_combo_allCells_byTrial_100X.pdf'),'-dpdf','-fillpage')
        [p_like0, table_like0, stats_like0] = anovan(likely_all,{off_id});
    end
end

maxloglike_change_diff = cell(ndiff,noff_all,nboot+1);
for iboot = 1:nboot+1
    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        for ioff = 1:noff_all
            for i =1:length(del)
                maxloglike_change_diff{idiff,ioff,iboot} = [maxloglike_change_diff{idiff,ioff,iboot} maxloglike_change{del(i),ioff,iboot}];
            end
        end
    end
end
save(fullfile(LG_base_out, 'maxloglike_change.mat'),'maxloglike_change_diff','maxloglike_change','loglike_all_fun')

figure;
dv = zeros(91,ndiff,noff_all,nboot+1);
for iboot = 1:nboot+1
    for idiff = 1:ndiff
        for ioff = 1:noff_all
            for itheta = 1:91
                dv(itheta,idiff,ioff,iboot) = length(find(maxloglike_change_diff{idiff,ioff,iboot}>(itheta-1)))./size(maxloglike_change_diff{idiff,ioff,iboot},2);
            end
        end
    end
end

dv_sort = sort(dv(:,:,:,2:end),4,'ascend');
dv_lb = dv(:,:,:,1)-dv_sort(:,:,:,5);
dv_ub = dv_sort(:,:,:,95)-dv(:,:,:,1);

figure;
subplot(2,2,1)
for ioff = 1:noff_all
    shadedErrorBar(1:91,dv(:,1,ioff,1),std(dv(:,1,ioff,2:end),[],4),{[col_mat(ioff,:) '-'], 'markerfacecolor',col_mat(ioff,:)},0)
    hold on
end
ylabel('FA rate')
xlabel('Decision variable (deg)')
ylim([0 1])

subplot(2,2,2)
for ioff = 1:noff_all
    shadedErrorBar(1:91,dv(:,2,ioff,1),std(dv(:,2,ioff,2:end),[],4),{[col_mat(ioff,:) '-'], 'markerfacecolor',col_mat(ioff,:)},0)
    hold on
end
ylabel('Hit rate')
xlabel('Decision variable (deg)')
ylim([0 1])
suptitle([mouse_list])
print(fullfile(LG_base_out, 'allCell_CDFP_decisionVariable.pdf'),'-dpdf','-fillpage')

maxloglike_change_diff_all = [];
maxloglike_change_diff_all_1 = [];
maxloglike_change_diff_all_2 = [];
off_id = [];
diff_id = [];
off_id_1 = [];
diff_id_1 = [];
off_id_2 = [];
diff_id_2 = [];
subplot(2,2,3)
for idiff = 1:ndiff
    for ioff = 1:noff_all
        errorbar(diffs(idiff),mean(maxloglike_change_diff{idiff,ioff,1},2), std(maxloglike_change_diff{idiff,ioff,1},[],2)./sqrt(size(maxloglike_change_diff{idiff,ioff,1},2)),['-o' col_mat(ioff,:)]);
        hold on
        sz = size(maxloglike_change_diff{idiff,ioff,1});
        maxloglike_change_diff_all = [maxloglike_change_diff_all maxloglike_change_diff{idiff,ioff,1}];
        off_id = [off_id ioff*ones(sz)];
        diff_id = [diff_id idiff*ones(sz)];
        if idiff == 1
            off_id_1 = [off_id_1 ioff*ones(sz)];
            diff_id_1 = [diff_id_1 idiff*ones(sz)];
            maxloglike_change_diff_all_1 = [maxloglike_change_diff_all_1 maxloglike_change_diff{idiff,ioff,1}];
        end
        if idiff == 2
            off_id_2 = [off_id_2 ioff*ones(sz)];
            diff_id_2 = [diff_id_2 idiff*ones(sz)];
            maxloglike_change_diff_all_2 = [maxloglike_change_diff_all_2 maxloglike_change_diff{idiff,ioff,1}];
        end
    end
end
xlabel('Stim - Adapter')
ylabel('Max likely ori - Adapter')
xlim([-10 100])
ylim([-10 100])
suptitle([mouse_list '- Max Log Likelihood dist from Adapter - by trial- scaled sub- 100X'])
print(fullfile(LG_base_out, 'maxloglike_distAdapt_byTrial_scaledsub_100X.pdf'),'-dpdf','-fillpage')
% [p_max, table_max, stats_max] = anovan(maxloglike_change_diff_all,{off_id, diff_id});
% [p_max_1, table_max_1, stats_max_1] = anovan(maxloglike_change_diff_all_1,{off_id_1});
% [p_max_2, table_max_2, stats_max_2] = anovan(maxloglike_change_diff_all_2,{off_id_2});
% [comp_max_1] = multcompare(stats_max_1);
% [comp_max_2] = multcompare(stats_max_2);



figure;
loglike_diff_all = [];
loglike_diff_all_1 = [];
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    for ioff = 1:noff_all
        loglike_diff{idiff,ioff} = [];
        for idel = 1:length(del)
            loglike_diff{idiff,ioff} = [loglike_diff{idiff,ioff} loglike_all_fun{del(idel),ioff,1}];
        end
        temp_likely = loglike_diff{idiff,ioff};
        errorbar(diffs(idiff),mean(temp_likely(end,:),2), std(temp_likely(end,:),[],2)./sqrt(size(temp_likely,2)),['-o' col_mat(ioff,:)]);
        hold on
        sz = size(loglike_diff{idiff,ioff});
        loglike_diff_all = [loglike_diff_all temp_likely(end,:)];
        if idiff == 1
            loglike_diff_all_1 = [loglike_diff_all_1 temp_likely(end,:)];
        end
    end
end
xlabel('Stim - Adapter')
ylabel('Likelyhood = 0')
xlim([-10 100])
ylim([-18000 0])
suptitle([mouse_list '- Log Likelihood is 0 - by trial- scaled sub- 100X'])
print(fullfile(LG_base_out, 'loglike_is0_distAdapt_byTrial_scaledsub_100X.pdf'),'-dpdf','-fillpage')
% [p_like0, table_like0, stats_like0] = anovan(loglike_diff_all,{off_id, diff_id});
% [p_like0_1, table_like0_1, stats_like0_1] = anovan(loglike_diff_all_1,{off_id_1});
% [comp_like0_1] = multcompare(stats_like0_1);
% %save(fullfile(LG_base_out, 'loglike_params.mat'),'loglike_all_neurons','loglike_all_sum','loglike_all_fact')

%% max log likelihood by trials by experiment
deltas = 22.5:22.5:180;
delta_diff = 180-deltas;
delta_diff(find(delta_diff>90)) = 180- delta_diff(find(delta_diff>90));
diffs = unique(delta_diff);
ndiff = length(diffs);

ppResp_trials = cell(nDelta,noff_all,nexp);
loglike_all_neurons = cell(nDelta,noff_all,nexp);
loglike_all_sum = cell(nDelta,noff_all,nexp);
loglike_all_fact = cell(nDelta,noff_all,nexp);
loglike_all_fun = cell(nDelta,noff_all,nexp);
maxloglike_all = cell(nDelta,noff_all,nexp);
maxloglike_change = cell(nDelta,noff_all,nexp);
maxloglike_change_dv = cell(nDelta,noff_all,nexp);
maxloglike_change_diff = cell(ndiff,noff_all,nexp);
maxloglike_change_diff_dv = cell(ndiff,noff_all,nexp);
max_ind = cell(1,nexp);
good_ind = cell(1,nexp);
sz_fit = size(fit_all,1);
cellN = zeros(1,nexp);
start = 1;
for iexp = 1:nexp
    ppTemp = ppResp_all{iexp};
    nc = size(ppTemp{1,1},1);
    good_ind{1,iexp} = find(theta_90_all(start:nc+start-1,:)<=22.5);
    cellN(1,iexp) = length(good_ind{1,iexp});
    [max_fit, max_ind{1,iexp}] = max(fit_all(:,start:nc+start-1,3),[],1);
    fit_all_thresh =  bsxfun(@rdivide,fit_all(:,start:nc+start-1,3),max_fit)*100;
    fit_all_thresh(find(fit_all_thresh<1)) = 1;
    fit_all_thresh(find(fit_all_thresh>170)) = 170;
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            ppT = ppTemp{ioff,idelta};
            nTrials = size(ppT,2);
            ppT_norm = bsxfun(@rdivide, ppT, max_fit')*100;
            ppT_norm(find(ppT_norm<1)) = 1;
            ppT_norm(find(ppT_norm>170)) = 170;
            ppResp_trials{idelta,ioff,iexp} = ppT_norm;

            loglike_all_neurons{idelta,ioff,iexp} = zeros(sz_fit,cellN(1,iexp),nTrials);
            loglike_all_sum{idelta,ioff,iexp} = zeros(sz_fit,cellN(1,iexp),nTrials);
            loglike_all_fact{idelta,ioff,iexp} = zeros(1,cellN(1,iexp),nTrials);
            for i = 1:nTrials
                for iCell = 1:cellN(1,iexp)
                    iC = good_ind{1,iexp}(iCell);
                    loglike_all_neurons{idelta,ioff,iexp}(:,iCell,i) = squeeze(bsxfun(@times,log(fit_all_thresh(:,iC)'),ppT_norm(iC,i)))';
                    loglike_all_sum{idelta,ioff,iexp}(:,iCell,i) = fit_all_thresh(:,iC);
                    loglike_all_fact{idelta,ioff,iexp}(1,iCell,i) = log(gamma(ppT_norm(iC,i)+1));
                end
            end
%             loglike_all_neurons{idelta,ioff,iexp} = zeros(sz_fit,cellN(1,iexp));
%             loglike_all_sum{idelta,ioff,iexp} = zeros(sz_fit,cellN(1,iexp));
%             loglike_all_fact{idelta,ioff,iexp} = zeros(1,cellN(1,iexp));
%             for i = 1:nTrials
%                 for iCell = 1:cellN(1,iexp)
%                     iC = good_ind{1,iexp}(iCell);
%                     loglike_all_neurons{idelta,ioff,iexp}(:,iCell) = squeeze(bsxfun(@times,log(fit_all_thresh(:,iC)'),mean(ppResp_trials{idelta,ioff,iexp}(iC,:),2)))';
%                     loglike_all_sum{idelta,ioff,iexp}(:,iCell) = fit_all_thresh(:,iC);
%                     loglike_all_fact{idelta,ioff,iexp}(1,iCell) = log(gamma(mean(ppResp_trials{idelta,ioff,iexp}(iC,:),2)+1));
%                 end
%             end
            loglike_all_fun{idelta,ioff,iexp} =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons{idelta,ioff,iexp},2), bsxfun(@plus,(nansum(loglike_all_sum{idelta,ioff,iexp},2)),nansum(loglike_all_fact{idelta,ioff,iexp},2))));
            [max_val maxloglike_all{idelta,ioff,iexp}] = max(loglike_all_fun{idelta,ioff,iexp},[],1);
            maxloglike_change{idelta,ioff,iexp} = maxloglike_all{idelta,ioff,iexp}-1;
            maxloglike_change{idelta,ioff,iexp}(find(maxloglike_all{idelta,ioff,iexp}>90)) = 180-maxloglike_change{idelta,ioff,iexp}(find(maxloglike_all{idelta,ioff,iexp}>90));
            maxloglike_change_dv{idelta,ioff,iexp}=zeros(1,91);
            for itheta = 0:90
                maxloglike_change_dv{idelta,ioff,iexp}(1,itheta+1) = length(find(maxloglike_change{idelta,ioff,iexp}>itheta))/nTrials;
            end
        end
    end
    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        for ioff = 1:noff_all
            for i =1:length(del)
                maxloglike_change_diff{idiff,ioff,iexp} = [maxloglike_change_diff{idiff,ioff,iexp} maxloglike_change{del(i),ioff,iexp}];
                for itheta = 0:90
                    maxloglike_change_diff_dv{idiff,ioff,iexp}(1,itheta+1) = length(find(maxloglike_change_diff{idiff,ioff,iexp}>itheta))/size(maxloglike_change_diff{idiff,ioff,iexp},2);
                end
            end
            maxloglike_change_diff{idiff,ioff,iexp} = mean(maxloglike_change_diff{idiff,ioff,iexp},2);
        end
    end
    start = start+nc;
end
    
avgFADV = zeros(91,nexp,noff);
avgCDDV = zeros(91,nexp,noff);
for iexp = 1:nexp
    for ioff = 1:noff_all
        avgFADV(:,iexp,ioff) = maxloglike_change_dv{8,ioff,iexp}';
        avgCDDV(:,iexp,ioff) = maxloglike_change_diff_dv{2,ioff,iexp}';
    end
end
avgReadout = zeros(ndiff,nexp,noff);
figure;
for iexp = 1:nexp
    subplot(3,4,iexp)
    for ioff = 1:noff_all
        for idiff = 1:ndiff
            avgReadout(idiff,iexp,ioff) = mean(maxloglike_change_diff{idiff,ioff,iexp},2);
            %avgReadout(idiff,iexp,ioff) = maxloglike_change_diff{idiff,ioff,iexp};
        end
        plot(diffs,avgReadout(:,iexp,ioff),'-o')
        hold on
    end
    xlim([0 90])
    ylim([0 90])
    title(num2str(cellN(1,iexp)))
end
figure;
good_exp = find(cellN>15);
subplot(2,2,1)
for ioff = 1:noff_all
    shadedErrorBar(0:90,mean(avgCDDV(:,good_exp,ioff),2),std(avgCDDV(:,good_exp,ioff),[],2)./sqrt(length(good_exp)),{['-' col_mat(ioff,:)],'markerfacecolor',col_mat(ioff,:)},0);
    hold on
end
ylabel('Hit Rate')
xlabel('Decision variable (deg)')
subplot(2,2,2)
for ioff = 1:noff_all
    shadedErrorBar(0:90,mean(avgFADV(:,good_exp,ioff),2),std(avgFADV(:,good_exp,ioff),[],2)./sqrt(length(good_exp)),{['-' col_mat(ioff,:)],'markerfacecolor',col_mat(ioff,:)},0);
    hold on
end
ylabel('FA Rate')
xlabel('Decision variable (deg)')
subplot(2,2,3)
for ioff = 1:noff_all
    errorbar(diffs, mean(avgReadout(:,good_exp,ioff),2),std(avgReadout(:,good_exp,ioff),[],2)./sqrt(length(good_exp)),'-o');
    hold on
end
xlim([0 100])
ylim([0 100])
xlabel('Stimulus (deg)')
ylabel('Readout (deg)')
suptitle([reshape(flipud(rot90(mouse_mat(good_exp,:))),[1 size(good_exp,2)*4]) '- Hit and FA rate by experiment by trial'])
print(fullfile(LG_base_out, 'FACD_dist_byExpt_byTrial.pdf'),'-dpdf','-fillpage')

%% max log likelihood by trials - exponential distribution
OSI_cells = good_ind_theta;
cellN = length(OSI_cells);
sz_fit = size(fit_all,1);
n_ind = zeros(nexp, nDelta,noff_all);
for iexp = 1:nexp
    ppTemp = ppResp_all{iexp};
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            n_ind(iexp,idelta,ioff) = size(ppTemp{ioff, idelta},2);
        end
    end
end
min_trials = squeeze(min(n_ind,[],1));

[max_fit, max_ind] = max(fit_all(:,:,3),[],1);
fit_all_thresh =  bsxfun(@rdivide,fit_all(:,:,3),max_fit)*100;
fit_all_thresh(find(fit_all_thresh<1)) = 1;
fit_all_thresh(find(fit_all_thresh>170)) = 170;

ppResp_trials = cell(noff_all,nDelta);
start = 1;
for iexp = 1:nexp
    ppTemp = ppResp_all{iexp};
    nc = size(ppTemp{1,1},1);
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            ppT = ppTemp{ioff,idelta};
            ppT_norm = bsxfun(@rdivide, ppT(:,1:min_trials(idelta,ioff)), max_fit(1,start:nc+start-1)')*100;
            ppT_norm(find(ppT_norm<1)) = 1;
            ppT_norm(find(ppT_norm>170)) = 170;
            ppResp_trials{ioff,idelta}(start:nc+start-1,:) = ppT_norm;
        end
    end
    start = start+nc;
end

cellN = length(good_ind_theta);
nCells = size(fit_all_thresh,2);
nboot = 0;
sz_fit = size(fit_all,1);
loglike_all_neurons = cell(nDelta,noff_all,nboot+1);
loglike_all_sum = cell(nDelta,noff_all,nboot+1);
loglike_all_fun = cell(nDelta,noff_all,nboot+1);
maxloglike_all = cell(nDelta,noff_all,nboot+1);
maxloglike_change = cell(nDelta,noff_all,nboot+1);
maxloglike_change_diff = cell(ndiff,noff_all,nboot+1);
for iboot = 1:nboot+1
    if iboot>1
        ind_cells_temp = good_ind_theta(randsample(cellN, cellN, 1));
    else
        ind_cells_temp = good_ind_theta;
    end
    fprintf('.')
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            loglike_all_neurons{idelta,ioff,iboot} = zeros(sz_fit,cellN,min_trials(idelta,ioff));
            loglike_all_fit{idelta,ioff,iboot} = zeros(sz_fit,cellN);
            for iCell = 1:cellN
                iC = ind_cells_temp(iCell);
                for i = 1:min_trials(idelta,ioff)
                    loglike_all_neurons{idelta,ioff,iboot}(:,iC,i) = squeeze(bsxfun(@times,sqrt(fit_all_thresh(:,iC))',ppResp_trials{ioff,idelta}(iC,i)))';
                end
                loglike_all_fit{idelta,ioff,iboot}(:,iC) = fit_all_thresh(:,iC);
            end
            loglike_all_fun{idelta,ioff,iboot} =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons{idelta,ioff,iboot}(:,ind_cells_temp,:),2), nansum(loglike_all_fit{idelta,ioff,iboot}(:,ind_cells_temp,:),2)));
            [max_val maxloglike_all{idelta,ioff,iboot}] = max(loglike_all_fun{idelta,ioff,iboot},[],1);
            maxloglike_change{idelta,ioff,iboot} = maxloglike_all{idelta,ioff,iboot};
            maxloglike_change{idelta,ioff,iboot}(find(maxloglike_all{idelta,ioff,iboot}>90)) = 179-maxloglike_change{idelta,ioff,iboot}(find(maxloglike_all{idelta,ioff,iboot}>90));
        end
    end

    if iboot == 1
        figure;
        for idelta = 1:nDelta
            for ioff = 3%1:noff_all
                subplot(3,3,idelta)
                shadedErrorBar(0:180, mean(loglike_all_fun{idelta,ioff,iboot},2), std(loglike_all_fun{idelta,ioff,iboot},[],2)./sqrt(size(loglike_all_fun{idelta,ioff,iboot},2)),{['-c'],'markerfacecolor','c'});
                hold on
            end
            title([num2str(deltas(idelta)) 'deg'])
            xlabel('Orientation')
            ylabel('Log likelihood')
            xlim([-10 190])
            %ylim([-3000 0])
        end
        suptitle([mouse_list '- Log Likelihood Function- all fits- by trial- 100X- sqrt'])
        print(fullfile(LG_base_out, 'sqrt_loglike_byDelta_byInt_summary_allFits_byTrial.pdf'),'-dpdf','-fillpage')

        %likelihood is 0
        figure;
        likely_all = [];
        off_id = [];
        for ioff = 1:noff_all
            for idelta = 1:nDelta
                temp_likely = loglike_all_fun{idelta,ioff,iboot};
                errorbar(deltas(idelta), mean(temp_likely(end,:),2),std(temp_likely(end,:),[],2)./sqrt(size(temp_likely,2)),['o' col_mat(ioff,:)]);
                hold on
                if idelta == nDelta
                    errorbar(0, mean(temp_likely(end,:),2),std(temp_likely(end,:),[],2)./sqrt(size(temp_likely,2)),['o' col_mat(ioff,:)]);
                end
                if idelta == 1
                    likely_all = [likely_all temp_likely(end,:)];
                    off_id = [off_id ioff*ones(size(temp_likely(1,:)))];
                end
            end
            xlim([-22 202])
            xlabel('Orientation')
            ylabel('Log likelihood of 0 deg')
        end
        suptitle([mouse_list '- Log Likelihood of 0 deg stimulus- all cells- by trial- 100X- sqrt'])
        print(fullfile(LG_base_out, 'sqrt_loglike_is0_combo_allCells_byTrial_100X.pdf'),'-dpdf','-fillpage')
        [p_like0, table_like0, stats_like0] = anovan(likely_all,{off_id});
    end

    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        for ioff = 1:noff_all
            for i =1:length(del)
                maxloglike_change_diff{idiff,ioff,iboot} = [maxloglike_change_diff{idiff,ioff,iboot} maxloglike_change{del(i),ioff,iboot}];
            end
        end
    end
end
save(fullfile(LG_base_out, 'maxloglike_change_sqrt.mat'),'maxloglike_change_diff','maxloglike_change','loglike_all_fun')

figure;
maxloglike_change_diff_all = [];
maxloglike_change_diff_all_1 = [];
maxloglike_change_diff_all_2 = [];
off_id = [];
diff_id = [];
off_id_1 = [];
diff_id_1 = [];
off_id_2 = [];
diff_id_2 = [];
for idiff = 1:ndiff
    for ioff = 3 %1:noff_all
        errorbar(diffs(idiff),mean(maxloglike_change_diff{idiff,ioff,1},2), std(maxloglike_change_diff{idiff,ioff,1},[],2)./sqrt(size(maxloglike_change_diff{idiff,ioff,1},2)), 'ok')%['-o' col_mat(ioff,:)]);
        hold on
        sz = size(maxloglike_change_diff{idiff,ioff,1});
        maxloglike_change_diff_all = [maxloglike_change_diff_all maxloglike_change_diff{idiff,ioff,1}];
        off_id = [off_id ioff*ones(sz)];
        diff_id = [diff_id idiff*ones(sz)];
        if idiff == 1
            off_id_1 = [off_id_1 ioff*ones(sz)];
            diff_id_1 = [diff_id_1 idiff*ones(sz)];
            maxloglike_change_diff_all_1 = [maxloglike_change_diff_all_1 maxloglike_change_diff{idiff,ioff,1}];
        end
        if idiff == 2
            off_id_2 = [off_id_2 ioff*ones(sz)];
            diff_id_2 = [diff_id_2 idiff*ones(sz)];
            maxloglike_change_diff_all_2 = [maxloglike_change_diff_all_2 maxloglike_change_diff{idiff,ioff,1}];
        end
    end
end
xlabel('Stim - Adapter')
ylabel('Max likely ori - Adapter')
xlim([-10 100])
ylim([-10 100])
axis square
suptitle([mouse_list '- Max Log Likelihood dist from Adapter - by trial- scaled sub- 100X - sqrt'])
print(fullfile(LG_base_out, 'sqrt_maxloglike_distAdapt_byTrial_scaledsub_100X.pdf'),'-dpdf','-fillpage')
[p_max, table_max, stats_max] = anovan(maxloglike_change_diff_all,{off_id, diff_id});
[p_max_1, table_max_1, stats_max_1] = anovan(maxloglike_change_diff_all_1,{off_id_1});
[p_max_2, table_max_2, stats_max_2] = anovan(maxloglike_change_diff_all_2,{off_id_2});
[comp_max_1] = multcompare(stats_max_1);
[comp_max_2] = multcompare(stats_max_2);

figure;
loglike_diff_all = [];
loglike_diff_all_1 = [];
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    for ioff = 1:noff_all
        loglike_diff{idiff,ioff} = [];
        for idel = 1:length(del)
            loglike_diff{idiff,ioff} = [loglike_diff{idiff,ioff} loglike_all_fun{del(idel),ioff,1}];
        end
        temp_likely = loglike_diff{idiff,ioff};
        errorbar(diffs(idiff),mean(temp_likely(end,:),2), std(temp_likely(end,:),[],2)./sqrt(size(temp_likely,2)),['-o' col_mat(ioff,:)]);
        hold on
        sz = size(loglike_diff{idiff,ioff});
        loglike_diff_all = [loglike_diff_all temp_likely(end,:)];
        if idiff == 1
            loglike_diff_all_1 = [loglike_diff_all_1 temp_likely(end,:)];
        end
    end
end
xlabel('Stim - Adapter')
ylabel('Likelyhood = 0')
xlim([-10 100])
ylim([-10000 0])
suptitle([mouse_list '- Log Likelihood is 0 - by trial- scaled sub- 100X - sqrt'])
print(fullfile(LG_base_out, 'sqrt_loglike_is0_distAdapt_byTrial_scaledsub_100X.pdf'),'-dpdf','-fillpage')
[p_like0, table_like0, stats_like0] = anovan(loglike_diff_all,{off_id, diff_id});
[p_like0_1, table_like0_1, stats_like0_1] = anovan(loglike_diff_all_1,{off_id_1});
[comp_like0_1] = multcompare(stats_like0_1);
%save(fullfile(LG_base_out, 'loglike_params.mat'),'loglike_all_neurons','loglike_all_sum','loglike_all_fact')

%% ROC
h = zeros(nDelta,1);
p = zeros(nDelta,1);
for idelta = 1:nDelta
    [h(idelta,:),p(idelta,:)] = ttest(roc_resp_all(:,1,idelta),roc_resp_all(:,2,idelta));
end
h_roc = zeros(nDelta,nDelta);
p_roc = zeros(nDelta,nDelta);
for i = 1:nDelta
    for idel = 1:nDelta
        ind_cells = intersect(good_ind_theta, find(max_dir_all == idel));
        [h_roc(idel,i), p_roc(idel,i)] = ttest(roc_resp_all(ind_cells,1,i),roc_resp_all(ind_cells,2,i));
    end
end

h_roc_diff = zeros(ndiff,ndiff);
p_roc_diff = zeros(ndiff,ndiff);
for i = 1:ndiff
    idiff = find(delta_diff == diffs(i)); 
    for idel = 1:ndiff
        dels = find(delta_diff == diffs(idel));
        ind_cells = [];
        for ii = 1:length(dels)
            ind_cells = [ind_cells; intersect(good_ind_theta, find(max_dir_all == dels(ii)))];
        end
        [h_roc_diff(idel,i), p_roc_diff(idel,i)] = ttest(mean(roc_resp_all(ind_cells,1,idiff),3),mean(roc_resp_all(ind_cells,2,idiff),3));
    end
end

figure
avg_roc = zeros(2,2,ndiff);
roc_all = [];
idiff_all = [];
ioff_all = [];
for idiff = 1:ndiff
    del = find(delta_diff ==  diffs(idiff));
    subplot(2,3,idiff)
    scatter(mean(roc_resp_all(good_ind_theta,1,del),3), mean(roc_resp_all(good_ind_theta,2,del),3), 'ob')
    roc_all = [roc_all ([mean(roc_resp_all(good_ind_theta,1,del),3); mean(roc_resp_all(good_ind_theta,2,del),3)])];
    avg_roc(:,1,idiff) = [mean(mean(roc_resp_all(good_ind_theta,1,del),3),1); mean(mean(roc_resp_all(good_ind_theta,2,del),3),1)];
    avg_roc(:,2,idiff) = [std(mean(roc_resp_all(good_ind_theta,1,del),3),[],1)./sqrt(length(good_ind_theta)); std(mean(roc_resp_all(good_ind_theta,2,del),3),[],1)./sqrt(length(good_ind_theta))];
    [htemp ptemp] = ttest(mean(roc_resp_all(good_ind_theta,1,del),3), mean(roc_resp_all(good_ind_theta,2,del),3));
    axis square
    xlim([0 1])
    ylim([0 1])
    refline(1,0)
    vline(0.5)
    hline(0.5)
    xlabel([num2str(chop(off_all(1,:)*(1000/frameRateHz),2)) ' ms ISI'])
    ylabel([num2str(chop(off_all(2,:)*(1000/frameRateHz),2)) ' ms ISI'])
    title(['0 vs ' num2str(diffs(idiff)) '- p = ' num2str(chop(ptemp,2))])
end
subplot(2,3,idiff+1)
errorbar(repmat(diffs, [2 1])', squeeze(avg_roc(:,1,:))', squeeze(avg_roc(:,2,:))', '-o')
ylim([0 1])
xlim([-10 100])
axis square
ylabel('Average auROC')
xlabel('Diff of Stim from Adaptor (deg)')
suptitle([mouse_list '- ROC 0 (short ISI) vs Diff from Adapter']) 
print(fullfile(LG_base_out, 'rocv180_allCells_summary.pdf'),'-dpdf','-fillpage')
[p_roc_all, table_roc_all, stats_roc_all] = anova2(roc_all, length(good_ind_theta));
comp_roc_all = multcompare(stats_roc_all);

figure
avg_roc_abs = zeros(2,2,ndiff);
for idiff = 1:ndiff
    del = find(delta_diff ==  diffs(idiff));
    subplot(2,3,idiff)
    scatter(abs(mean(roc_resp_all(good_ind_theta,1,del),3)-0.5), abs(mean(roc_resp_all(good_ind_theta,2,del),3)-0.5), 'ob')
    avg_roc_abs(:,1,idiff) = [mean(abs(mean(roc_resp_all(good_ind_theta,1,del),3)-0.5),1); mean(abs(mean(roc_resp_all(good_ind_theta,2,del),3)-0.5),1)];
    avg_roc_abs(:,2,idiff) = [std(abs(mean(roc_resp_all(good_ind_theta,1,del),3)-0.5),[],1)./sqrt(length(good_ind_theta)); std(abs(mean(roc_resp_all(good_ind_theta,2,del),3)-0.5),[],1)./sqrt(length(good_ind_theta))];
    [htemp ptemp] = ttest(abs(mean(roc_resp_all(good_ind_theta,1,del),3)-0.5), abs(mean(roc_resp_all(good_ind_theta,2,del),3)-0.5));
    axis square
    xlim([0 .5])
    ylim([0 .5])
    refline(1,0)
    vline(0.5)
    hline(0.5)
    xlabel([num2str(chop(off_all(1,:)*(1000/frameRateHz),2)) ' ms ISI'])
    ylabel([num2str(chop(off_all(2,:)*(1000/frameRateHz),2)) ' ms ISI'])
    title(['0 vs ' num2str(diffs(idiff)) '- p = ' num2str(chop(ptemp,2))])
end
subplot(2,3,idiff+1)
errorbar(repmat(diffs, [2 1])', squeeze(avg_roc_abs(:,1,:))', squeeze(avg_roc_abs(:,2,:))', '-o')
ylim([0 0.5])
xlim([-10 100])
axis square
ylabel('Average auROC')
xlabel('Diff of Stim from Adaptor (deg)')
suptitle([mouse_list '- ROC 0 (short ISI) vs Diff from Adapter- Absolute value']) 
print(fullfile(LG_base_out, 'rocv180_allCells_summary_abs.pdf'),'-dpdf','-fillpage')

%bin ROC by ori pref
figure;
bins = 0:15:75;
nbin = length(bins);
avg_roc_bin = zeros(2,2,ndiff,nbin);
[n,n2] = subplotn(nbin);
for ibin = 1:nbin
    ind = intersect(good_ind_theta,find(pref_ori_all_diff>bins(ibin)));
    for idiff = 1:ndiff
        del = find(delta_diff ==  diffs(idiff));
        avg_roc_bin(:,1,idiff,ibin) = [mean(mean(roc_resp_all(ind,1,del),3),1); mean(mean(roc_resp_all(ind,2,del),3),1)];
        avg_roc_bin(:,2,idiff,ibin) = [std(mean(roc_resp_all(ind,1,del),3),[],1)./sqrt(length(ind)); std(mean(roc_resp_all(ind,2,del),3),[],1)./sqrt(length(ind))];
    end
    subplot(n,n2,ibin)
    errorbar([diffs' diffs'], squeeze(avg_roc_bin(:,1,:,ibin))',squeeze(avg_roc_bin(:,2,:,ibin))','-o');
    title([num2str(bins(ibin)) '-90deg pref - n = ' num2str(length(ind))])
    ylim([0.5 1])
    xlabel('Test-Adapter')
    ylabel('auROC')
    legend(strvcat('250', '750'),'Location','Northwest')
end
suptitle([mouse_list '- ROC 0 (short ISI) vs Diff from Adapter- Binned Cell Groups']) 
print(fullfile(LG_base_out, 'roc_binnedByPref_summary.pdf'),'-dpdf','-fillpage')

%bin with pos and neg weights
figure;
bins = 0:15:91;
nbin = length(bins);
roc_val = nan(nDelta,noff,nexp,nbin,nbin);
roc_val_diff = nan(ndiff,noff,nexp,nbin,nbin);
for ibin = 1:nbin
    for ibin2 = 1:nbin
        offset = 0;
        t = 0;
        for iexp = 1:nexp
            ppResp = ppResp_all{iexp};
            n = size(ppResp{1,1},1);
            ind1 = find(theta_90_all(offset+1:offset+n,control_ind)<=22.5);
            ind2 = intersect(ind1,find(pref_ori_all_diff(offset+1:offset+n,control_ind)>bins(ibin)));
            ind3 = intersect(ind1,find(pref_ori_all_diff(offset+1:offset+n,control_ind)<bins(ibin2)));
            ind = intersect(ind2,ind3);
            if length(ind) == 0 & length(ind1)>9 & length([ind2; ind3])>0
                for iori = 1:nDelta
                    for ioff = 1:noff
                        if length(ind2)>0
                            noise_90 = mean(ppResp{1,8}(ind2,:),1);
                            sig_90 = mean(ppResp{ioff,iori}(ind2,:),1);
                        else
                            noise_90 = zeros(1,size(ppResp{1,8},2));
                            sig_90 = zeros(1,size(ppResp{ioff,iori},2));
                        end
                        if length(ind3)>0
                            noise_0 = mean(ppResp{1,8}(ind3,:),1);
                            sig_0 = mean(ppResp{ioff,iori}(ind3,:),1);
                        else
                            noise_0 = zeros(1,size(ppResp{1,8},2));
                            sig_0 = zeros(1,size(ppResp{ioff,iori},2));
                        end
                        roc_val(iori,ioff,iexp,ibin,ibin2) = roc_gh(noise_90-noise_0, sig_90-sig_0);
                    end
                end
                t = t+1;
            else
                roc_val(:,:,iexp,ibin,ibin2) = nan(nDelta,noff);
            end
            offset = n+offset;
        end
        for idiff = 1:ndiff
            del = find(delta_diff ==  diffs(idiff));
            roc_val_diff(idiff,:,:,ibin,ibin2) = nanmean(roc_val(del,:,:,ibin,ibin2),1);
        end
        subplot(nbin,nbin,ibin2+(nbin*(ibin-1)))
        for ioff = 1:noff
            errorbar(diffs, nanmean(roc_val_diff(:,ioff,:,ibin,ibin2),3), nanstd(roc_val_diff(:,ioff,:,ibin,ibin2),[],3)./sqrt(nexp), '-o')
            hold on
        end
        ylim([0 1])
        xlim([-10 100])
        if ibin == 1
            title(['Neg < ' num2str(bins(ibin2))])
        end
        if ibin2 == 1
            ylabel(['Pos > ' num2str(bins(ibin))])
        end
        text(0,0.1,['nexp = ' num2str(t)])
        %ylabel('auROC')
        %xlabel('Diff of Stim from Adaptor (deg)')
    end
end
suptitle('Blue- 250 ms; Red- 750 ms')
print(fullfile(LG_base_out, 'roc_Pos&NegWeight_summary.pdf'),'-dpdf','-fillpage')
%all pos 90 vs >30 pos 90
[h p] = ttest(roc_val_diff(end,2,:,1,1),roc_val_diff(end,2,:,3,1));
%all pos 90 vs >60 pos 90
[h p] = ttest(roc_val_diff(end,2,:,1,1),roc_val_diff(end,2,:,5,1));
%>30 pos 90vs >60 pos 90
[h p] = ttest(roc_val_diff(end,2,:,3,1),roc_val_diff(end,2,:,5,1));
%250 vs 750 for >30 pos and <15 neg for 22.5
[h p] = ttest(roc_val_diff(2,1,:,3,2),roc_val_diff(2,2,:,3,2));
%250 vs 750 for >30 pos and <15 neg for 0
[h p] = ttest(roc_val_diff(1,1,:,3,2),roc_val_diff(1,2,:,3,2));
%250 vs 750 for >30 pos and <30 neg for 22.5
[h p] = ttest(roc_val_diff(2,1,:,3,3),roc_val_diff(2,2,:,3,3));
%250 vs 750 for >30 pos and <30 neg for 0
[h p] = ttest(roc_val_diff(1,1,:,3,3),roc_val_diff(1,2,:,3,3));
%250 vs 750 for >60 pos and <15 neg for 22.5
[h p] = ttest(roc_val_diff(2,1,:,5,2),roc_val_diff(2,2,:,5,2));
%250 vs 750 for >60 pos and <15 neg for 0
[h p] = ttest(roc_val_diff(1,1,:,5,2),roc_val_diff(1,2,:,5,2));
%250 vs 750 for >60 pos and <30 neg for 22.5
[h p] = ttest(roc_val_diff(2,1,:,5,3),roc_val_diff(2,2,:,5,3));
%250 vs 750 for >60 pos and <30 neg for 0
[h p] = ttest(roc_val_diff(1,1,:,5,3),roc_val_diff(1,2,:,5,3));
%250 vs 750 for >60 pos for 22.5
[h p] = ttest(roc_val_diff(2,1,:,5,1),roc_val_diff(2,2,:,5,1));
%250 vs 750 for >60 pos for 0
[h p] = ttest(roc_val_diff(1,1,:,5,1),roc_val_diff(1,2,:,5,1));
%250 vs 750 for >30 pos for 22.5
[h p] = ttest(roc_val_diff(2,1,:,3,1),roc_val_diff(2,2,:,3,1));
%250 vs 750 for >30 pos for 0
[h p] = ttest(roc_val_diff(1,1,:,3,1),roc_val_diff(1,2,:,3,1));

%250 vs 750 for >60 pos
ind_nan = find(~isnan(roc_val_diff(1,1,:,5,1)));
[p,table] = anova2([squeeze(roc_val_diff(1:3,1,ind_nan,5,1))'; squeeze(roc_val_diff(1:3,2,ind_nan,5,1))'],length(ind_nan))
%250 vs 750 for >30 pos
ind_nan = find(~isnan(roc_val_diff(1,1,:,3,1)));
[p,table] = anova2([squeeze(roc_val_diff(1:3,1,ind_nan,3,1))'; squeeze(roc_val_diff(1:3,2,ind_nan,3,1))'],length(ind_nan))

%250 vs 750 for >30 pos
ind_nan = find(~isnan(roc_val_diff(1,1,:,3,1)));
[p1,table] = anova2([squeeze(roc_val_diff(2:5,1,ind_nan,3,1))'; squeeze(roc_val_diff(2:5,2,ind_nan,3,1))'],length(ind_nan));
%250 vs 750 for >60 pos
ind_nan = find(~isnan(roc_val_diff(1,1,:,5,1)));
[p2,table] = anova2([squeeze(roc_val_diff(2:5,1,ind_nan,5,1))'; squeeze(roc_val_diff(2:5,2,ind_nan,5,1))'],length(ind_nan));
%250 vs 750 for >30 pos & <15 neg
ind_nan = find(~isnan(roc_val_diff(1,1,:,3,2)));
[p3,table] = anova2([squeeze(roc_val_diff(2:5,1,ind_nan,3,2))'; squeeze(roc_val_diff(2:5,2,ind_nan,3,2))'],length(ind_nan));
%250 vs 750 for >60 pos & <15 neg
ind_nan = find(~isnan(roc_val_diff(1,1,:,5,2)));
[p4,table] = anova2([squeeze(roc_val_diff(2:5,1,ind_nan,5,2))'; squeeze(roc_val_diff(2:5,2,ind_nan,5,2))'],length(ind_nan));

%histograms for all experiments
offset = 0;
for iexp = 1:nexp
    figure;
    ppResp = ppResp_all{iexp};
    n = size(ppResp{1,1},1);
    ind = find(theta_90_all(offset+1:offset+n,:)<22.5);
    noise = mean(ppResp{1,8}(ind,:),1);
    sig_250 = mean(ppResp{1,1}(ind,:),1);
    sig_750 = mean(ppResp{2,1}(ind,:),1);
    minx = min([noise sig_250 sig_750],[],2);
    maxx = max([noise sig_250 sig_750],[],2);
    edges = minx:(maxx-minx)/10:maxx+(maxx-minx)/10;
    subplot(3,1,1)
    hist(noise,edges)
    xlim([(minx-((maxx-minx)/10)) (maxx+((maxx-minx)/10))])
    ylim([0 4])
    title(['Expt ' num2str(iexp)])
    subplot(3,1,2)
    hist(sig_250,edges)
    roc_val_250 = roc_gh(noise,sig_250);
    title([' 250 ms ISI- AuROC = ' num2str(chop(roc_val_250,2))])
    xlim([(minx-((maxx-minx)/10)) (maxx+((maxx-minx)/10))])
    ylim([0 4])
    subplot(3,1,3)
    hist(sig_750,edges)
    roc_val_750 = roc_gh(noise,sig_750);
    title([' 750 ms ISI- AuROC = ' num2str(chop(roc_val_750,2))])
    xlim([(minx-((maxx-minx)/10)) (maxx+((maxx-minx)/10))]) 
    ylim([0 4])
    offset = n+offset;
end
print(fullfile(LG_base_out, 'avgRespHist_expt10.pdf'),'-dpdf','-fillpage')
    
%average ROC for good cells in each experiment
offset = 0;
roc_val = nan(nDelta,noff,nexp);
n_15 = zeros(1,nexp);
for iexp = 1:nexp
    ppResp = ppResp_all{iexp};
    n = size(ppResp{1,1},1);
    ind = find(theta_90_all(offset+1:offset+n,:)<22.5);
    %if length(ind)>15
        for iori = 1:nDelta
            for ioff = 1:noff
                roc_val(iori,ioff,iexp) = roc_gh(mean(ppResp{1,8}(ind,:),1), mean(ppResp{ioff,iori}(ind,:),1));
            end
        end
        n_15(1,iexp) =length(ind);
    %end
    offset = n+offset;
end

figure;
subplot(2,2,1)
for ioff = 1:noff
    errorbar(deltas, nanmean(roc_val(:,ioff,:),3), nanstd(roc_val(:,ioff,:),[],3)./sqrt(nexp),'-o')
    hold on
end
ylim([0.5 1])
xlim([-10 190])
ylabel('auROC')
xlabel('Orientation (deg)')
%title(['theta 90 < 22.5; >15 cells; n = ' num2str(length(find(n_15>0))) 'expts'])

roc_val_diff = zeros(ndiff,noff,nexp);
for idiff = 1:ndiff
    del = find(delta_diff ==  diffs(idiff));
    roc_val_diff(idiff,:,:) = mean(roc_val(del,:,:),1);
end
subplot(2,2,2)
for ioff = 1:noff
    errorbar(diffs, nanmean(roc_val_diff(:,ioff,:),3), nanstd(roc_val_diff(:,ioff,:),[],3)./sqrt(nexp), '-o')
    hold on
end
ylim([0.5 1])
xlim([-10 100])
ylabel('auROC')
xlabel('Diff of Stim from Adaptor (deg)')
print(fullfile(LG_base_out, 'rocv180_allExpts_summary.pdf'),'-dpdf','-fillpage')
roc_val_diff_nonan = roc_val_diff;
roc_val_diff_nonan(:,:,find(n_15==0)) = [];
[p_roc_byexp, table_roc_byexp, stats_roc_byexp] = anova2(reshape(permute(roc_val_diff_nonan,[1 3 2]), [ndiff noff*length(find(n_15>0))])', length(find(n_15>0)));
for idiff = 1:ndiff
    [p_roc_byexp(idiff), h_roc_byexp(idiff)] = ttest(squeeze(roc_val_diff(idiff,1,:)),squeeze(roc_val_diff(idiff,2,:)));
end
strtmp = [];
for i = 2:4
    strtmp = [strtmp num2str(chop(table_roc_byexp{i,6},3)) ' '];
end
title(strtmp)

%hit/FA rate by DV
%FA rate
offset = 0;
nbin = 101;
dv_FA_mat = zeros(nbin,noff_all,nexp);
dv_CD_mat = zeros(nbin,noff_all,nexp);
figure
for iexp = 1:nexp
    subplot(3,4,iexp)
    ppResp = ppResp_all{iexp};
    for ioff = noff_all:-1:1
        data_resp_off_zero = ppResp{ioff,8};
        data_resp_off_22 = [ppResp{ioff,1} ppResp{ioff,7}];
        [n, ntrials_zero] = size(data_resp_off_zero);
        [n, ntrials_22] = size(data_resp_off_22);
        ind = find(theta_90_all(offset+1:offset+n,:)<22.5);
        data_resp_off_zero_avg = mean(data_resp_off_zero(ind,:),1);
        data_resp_off_22_avg = mean(data_resp_off_22(ind,:),1);
        if ioff == 3
            min_data_zero = min(data_resp_off_zero_avg,[],2);
            max_data_zero = max(data_resp_off_zero_avg,[],2);
            min_data_22 = min(data_resp_off_22_avg,[],2);
            max_data_22 = max(data_resp_off_22_avg,[],2);
            data_incr_zero = (max_data_zero-min_data_zero)./nbin;
            data_incr_22 = (max_data_22-min_data_22)./nbin;
        end
        for ibin = 1:nbin
            dv_FA_mat(ibin,ioff,iexp) = length(find(data_resp_off_zero_avg>=min_data_zero+((ibin-1)*(data_incr_zero))))./ntrials_zero;
            dv_CD_mat(ibin,ioff,iexp) = length(find(data_resp_off_22_avg>=min_data_22+((ibin-1)*(data_incr_22))))./ntrials_22;
        end
        plot(0:0.01:1,dv_FA_mat(:,ioff,iexp))
        hold on
    end
    title(length(ind))
    offset = n+offset;
end

figure; 
gray_mat = grays(3);
for ioff = 1:noff_all
    subplot(2,2,1)
    shadedErrorBar(0:0.01:1,mean(dv_CD_mat(:,ioff,:),3),std(dv_CD_mat(:,ioff,:),[],3)./sqrt(nexp),{col_mat(ioff,:),'markerfacecolor',col_mat(ioff,:)},0);
    hold on
    ylabel('CD rate')
    xlabel('Norm decision variable')
    subplot(2,2,2)
    shadedErrorBar(0:0.01:1,mean(dv_FA_mat(:,ioff,:),3),std(dv_FA_mat(:,ioff,:),[],3)./sqrt(nexp),{col_mat(ioff,:),'markerfacecolor',col_mat(ioff,:)},0);
    hold on
    ylabel('FA rate')
    xlabel('Norm decision variable')
    legend({'250','750','control'})
end

% for ioff = 1:noff_all
%     subplot(2,2,3)
%     errorbar(0:0.01:1,mean(dv_CD_mat(:,ioff,find(n_15)),3),std(dv_CD_mat(:,ioff,find(n_15)),[],3)./sqrt(length(find(n_15))));
%     hold on
%     ylabel('CD rate')
%     xlabel('Norm decision variable')
%     subplot(2,2,4)
%     errorbar(0:0.01:1,mean(dv_FA_mat(:,ioff,find(n_15)),3),std(dv_FA_mat(:,ioff,find(n_15)),[],3)./sqrt(length(find(n_15))));
%     hold on
%     ylabel('FA rate')
%     xlabel('Norm decision variable')
%     legend({'250','750','control'})
% end

print(fullfile(LG_base_out, 'CD_FA_sum_allExpts_summary.pdf'),'-dpdf','-fillpage')

%% tuning figures- sub
delta_resp_norm_all = bsxfun(@rdivide, delta_resp_sub_all, max(delta_resp_sub_all(:,:,3),[],2));
figure;
for idel = 1:nDelta
    ind_cells = intersect(good_ind_theta,find(max_dir_all == idel));
    for ioff = 1:noff_all
        subplot(noff_all,nDelta,((ioff-1)*nDelta)+idel)
        errorbar(deltas, mean(delta_resp_norm_all(ind_cells,:,ioff),1), std(delta_resp_norm_all(ind_cells,:,ioff),[],1)./sqrt(length(ind_cells)), 'ok')
        if ioff == 1
            title([num2str(deltas(idel)) ' deg pref- ' num2str(length(ind_cells)) ' cells'])
        end
        ylim([-0.1 1.2])
    end
end
suptitle([mouse_list '- avg all cells- normalized- sub'])
print(fullfile(LG_base_out, 'tuning_byPref_byInt_sep_summary_sub.pdf'),'-dpdf','-fillpage')

figure;
c = [.6 .6 .6; .4 .4 .4; 0 0 0];
[n n2] = subplotn(nDelta);
for idel = 1:nDelta
    subplot(n,n2,idel)
    ind_cells = intersect(good_ind_theta,find(max_dir_all == idel));
    for ioff = 1:noff_all
        errorbar(deltas, mean(delta_resp_norm_all(ind_cells,:,ioff),1), std(delta_resp_norm_all(ind_cells,:,ioff),[],1)./sqrt(length(ind_cells)), 'o', 'Color', c(ioff,:));
        hold on
        if ioff == 1
            title([num2str(deltas(idel)) ' deg pref- ' num2str(length(ind_cells)) ' cells'])
        end
        ylim([-0.1 1.2])
    end
    xlabel('Stimulus Orientation (deg)')
end
suptitle([mouse_list '- avg all cells- normalized- sub'])
print(fullfile(LG_base_out, 'tuning_byPref_byInt_summary_sub.pdf'),'-dpdf','-fillpage')

figure;
c = [.6 .6 .6; .4 .4 .4; 0 0 0];
[n n2] = subplotn(nDelta);
for idelta = 1:nDelta
    for ioff = 1:noff_all
        for idel = 1:nDelta
            ind_cells = intersect(good_ind_theta, find(max_dir_all == idel));
            subplot(n,n2,idelta)
            errorbar(deltas(idel), mean(delta_resp_norm_all(ind_cells,idelta,ioff),1),std(delta_resp_norm_all(ind_cells,idelta,ioff),[],1)./sqrt(length(ind_cells)),'o', 'Color', c(ioff,:));
            hold on
        end
    end
    ylim([-0.1 1.2])
    xlabel('Cell preference group (deg)')
    title([num2str(deltas(idelta)) ' deg change'])
end
suptitle([mouse_list '- avg all cells- normalized- sub'])
print(fullfile(LG_base_out, 'popTuning_byInt_summary_sub.pdf'),'-dpdf','-fillpage')


%% Pref and OSI change - sub
deltas = 22.5:22.5:180;
delta_diff = 180-deltas;
delta_diff(find(delta_diff>90)) = 180- delta_diff(find(delta_diff>90));
diffs = unique(delta_diff);
ndiff = length(diffs);

figure;
subplot(2,2,1)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_sub_all(cell_ind,ioff)-OSI_sub_all(cell_ind,3),1), std(OSI_sub_all(cell_ind,ioff)-OSI_sub_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in OSI'])
ylim([-.3 .3])
xlim([-10 100])
hline(0)
subplot(2,2,2)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_sub_all(cell_ind,ioff)./OSI_sub_all(cell_ind,3),1), std(OSI_sub_all(cell_ind,ioff)./OSI_sub_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Ratio of OSI'])
ylim([0.5 1.5])
xlim([-10 100])
hline(1)
subplot(2,2,3)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(OSI_k_sub_all(cell_ind,ioff)-OSI_k_sub_all(cell_ind,3),1), std(OSI_k_sub_all(cell_ind,ioff)-OSI_k_sub_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in OSI-k'])
title('OSI-k')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
subplot(2,2,4)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(HWHM_sub_all(cell_ind,ioff)-HWHM_sub_all(cell_ind,3),1), std(HWHM_sub_all(cell_ind,ioff)-HWHM_sub_all(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of max from adaptor (deg)'])
ylabel(['Difference in HWHM'])
title('HWHM')
ylim([-.3 .3])
xlim([-10 100])
hline(0)
suptitle([mouse_list '- OSI by Interval- Sub'])
print(fullfile(LG_base_out, 'OSI_byInt_summary_sub.pdf'),'-dpdf','-fillpage')

figure;
pref_ori_all_diff = pref_ori_sub_all;
pref_ori_all_diff(find(pref_ori_all_diff>90)) = 180-pref_ori_all_diff(find(pref_ori_all_diff>90));
subplot(2,2,1)
g = jet(180);
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    plot(pref_ori_sub_all(iC, 2),pref_ori_sub_all(iC, 1), 'o', 'Color', g(ceil(pref_ori_sub_all(iC,3))+1,:))
    hold on
end
refline(1,0)
axis square
xlabel(['Pref Ori- ' num2str(chop(off_all(2).*(1000/frameRateHz),3)) ' ms ISI'])
ylabel(['Pref Ori- ' num2str(chop(off_all(1).*(1000/frameRateHz),3)) ' ms ISI'])
subplot(2,2,2)
g = jet(91);
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    plot(pref_ori_all_diff(iC, 2),pref_ori_all_diff(iC, 1), 'o', 'Color', g(ceil(pref_ori_all_diff(iC,3)+1),:))
    hold on
end
refline(1,0)
axis square
xlabel(['Diff Pref from Adapt Ori- ' num2str(chop(off_all(2).*(1000/frameRateHz),3)) ' ms ISI'])
ylabel(['Diff Pref from Adapt Ori- ' num2str(chop(off_all(1).*(1000/frameRateHz),3)) ' ms ISI'])
subplot(2,2,3)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    cell_ind = [];
    for i = 1:length(del)
        cell_ind = [cell_ind; intersect(good_ind_theta, find(max_dir_all == del(i)))];
    end
    for ioff = 1:noff
        errorbar(diffs(idiff), mean(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),1), std(pref_ori_all_diff(cell_ind,ioff)-pref_ori_all_diff(cell_ind,3),[],1)./sqrt(length(cell_ind)), ['o' col_mat(ioff)])
        hold on
    end
end
xlabel(['Diff of Max from Adaptor (deg)'])
ylabel(['Change in Pref'])
ylim([-40 40])
hline(0)
xlim([-10 100])
suptitle([mouse_list '- Change in Preference by Interval- Sub'])
print(fullfile(LG_base_out, 'prefDiff_byInt_summary_sub.pdf'),'-dpdf','-fillpage')

col_mat = strvcat('b','r','y');
figure;
for i = 1:noff_all
    for idiff = 1:ndiff
        del = find(delta_diff == diffs(idiff));
        subplot(2,2,1)
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_sub_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_sub_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        hold on
        subplot(2,2,2)
        errorbar(diffs(idiff), squeeze(mean(mean(delta_resp_norm_all(good_ind_theta,del,i),2),1)), squeeze(std(mean(delta_resp_norm_all(good_ind_theta,del,i),2),[],1))./sqrt(length(good_ind_theta)), ['o' col_mat(i,:)])
        hold on
    end
end
subplot(2,2,1)
title('Absolute')
ylabel('Mean dF/F')
xlabel('Diff of Stim from Adaptor (deg)')
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([0 .3])
subplot(2,2,2)
title('Normalized')
ylabel('Mean dF/F')
xlabel('Diff of Stim from Adaptor (deg)')
set(gca, 'Xtick', 0:30:90)
xlim([-10 100])
ylim([0 1])
subplot(2,2,3)
for idiff = 1:ndiff
    del = find(delta_diff == diffs(idiff));
    diff_resp_group= [];
    for i = 1:length(del)
        cell_ind = intersect(good_ind_theta, find(max_dir_all == del(i)));
        diff_resp_group = [diff_resp_group; squeeze(delta_resp_norm_all(cell_ind,del(i),1:2))];
    end
    for i = 1:noff
        errorbar(diffs(idiff), mean(diff_resp_group(:,i),1), std(diff_resp_group(:,i),[],1)./sqrt(size(diff_resp_group,1)),['o' col_mat(i,:)])
        hold on
    end
end
title('Normalized')
ylabel('dF/F at max ori')
xlabel('Diff of Max from Adaptor (deg)')
ylim([0 2])
xlim([-10 100])
hline(1)
suptitle([mouse_list '- Mean Resp by Interval- sub'])
print(fullfile(LG_base_out, 'resp_byInt_summary_sub.pdf'),'-dpdf','-fillpage')

%% max log likelihood- sub

nboot = 1000;
OSI_cells = good_ind_theta;
cellN = length(OSI_cells);

loglike_all_neurons = zeros(cellN,nDelta,nDelta,noff_all);
loglike_all_fact = zeros(cellN,1,nDelta,noff_all);
loglike_all_sum = zeros(cellN,nDelta,nDelta,noff_all);
loglike_all_fun = zeros(nDelta,nDelta,noff_all,nboot+1);
loglike_all_fun_nosub = zeros(nDelta,nDelta,noff_all,nboot+1);
loglike_all_fun_subfact = zeros(nDelta,nDelta,noff_all,nboot+1);
maxloglike_all = zeros(nDelta,noff_all,nboot+1);
maxloglike_all_subfact = zeros(nDelta,noff_all,nboot+1);
maxloglike_all_nosub = zeros(nDelta,noff_all,nboot+1);

delta_resp_all_thresh = delta_resp_sub_all*100;
delta_resp_all_thresh(find(delta_resp_sub_all<0)) = 0.001;
for iCell = 1:length(good_ind_theta)
    iC = good_ind_theta(iCell);
    for idelta = 1:nDelta
        for ioff =1:noff_all
            loglike_all_neurons(iC,:,idelta,ioff) = squeeze(log10(delta_resp_all_thresh(iC,:,noff_all)).*delta_resp_all_thresh(iC,idelta,ioff))';
            loglike_all_sum(iC,:,idelta,ioff) = squeeze(delta_resp_all_thresh(iC,:,noff_all))';
            loglike_all_fact(iC,1,idelta,ioff) = log10(gamma(delta_resp_all_thresh(iC,idelta,ioff)+1));
        end
    end
end
    
for iboot = 1:nboot+1
    if iboot>1
        ind_cells_temp = OSI_cells(randsample(cellN, cellN, 1));
    else
        ind_cells_temp = OSI_cells;
    end

    loglike_all_fun(:,:,:,iboot) =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons(ind_cells_temp,:,:,:),1), bsxfun(@plus,(nansum(loglike_all_sum(ind_cells_temp,:,:,:),1)),nansum(loglike_all_fact(ind_cells_temp,:,:,:),1))));
    loglike_all_fun_subfact(:,:,:,iboot) =  squeeze(bsxfun(@minus, nansum(loglike_all_neurons(ind_cells_temp,:,:,:),1), nansum(loglike_all_fact(ind_cells_temp,:,:,:),1)));
    loglike_all_fun_nosub(:,:,:,iboot) =  squeeze(nansum(loglike_all_neurons(ind_cells_temp,:,:,:),1));
    for idelta = 1:nDelta
        for ioff = 1:noff_all
            [max_val maxloglike_all(idelta,ioff,iboot)] = max(loglike_all_fun(:,idelta,ioff,iboot),[],1);
            [max_val maxloglike_all_subfact(idelta,ioff,iboot)] = max(loglike_all_fun_subfact(:,idelta,ioff,iboot),[],1);
            [max_val maxloglike_all_nosub(idelta,ioff,iboot)] = max(loglike_all_fun_nosub(:,idelta,ioff,iboot),[],1);
        end
    end
end


figure;
for idelta = 1:nDelta
    subplot(3,3,idelta)
    for ioff = 1:noff_all
        plot(deltas, loglike_all_fun(:,idelta,ioff,1))
        hold on
    end
    title([num2str(deltas(idelta)) 'deg'])
    xlabel('Orientation')
    ylabel('Log likelihood')
end
suptitle([mouse_list '- Log Likelihood Function- all cells- 100X- Sub'])
print(fullfile(LG_base_out, 'loglike_byDelta_byInt_summary_allCells_100X_sub.pdf'),'-dpdf','-fillpage')



figure;
for ioff = 1:noff_all
    subplot(1,3,1)
    temp_likely = squeeze(loglike_all_fun_nosub(:,:,ioff,1));
    temp_sort = squeeze(sort(loglike_all_fun_nosub(:,:,ioff,:),4));
    errorbar([0 deltas], [temp_likely(end,end) temp_likely(end,:)],[squeeze(std(loglike_all_fun_nosub(end,end,ioff,2:end),[],4)) squeeze(std(loglike_all_fun_nosub(end,:,ioff,2:end),[],4))]);
    hold on
    if ioff == noff_all
        title('No sub')
    end
    subplot(1,3,2)
    temp_likely = squeeze(loglike_all_fun(:,:,ioff,1));
    temp_sort = squeeze(sort(loglike_all_fun(:,:,ioff,:),4));
    errorbar([0 deltas], [temp_likely(end,end) temp_likely(end,:)],[squeeze(std(loglike_all_fun(end,end,ioff,2:end),[],4)) squeeze(std(loglike_all_fun(end,:,ioff,2:end),[],4))]);
    
    hold on
    if ioff == noff_all
        title('Sub Fact and Sum')
    end
    
    subplot(1,3,3)
    temp_likely = squeeze(loglike_all_fun_subfact(:,:,ioff,1));
    temp_sort = squeeze(sort(loglike_all_fun_subfact(:,:,ioff,:),4));
    errorbar([0 deltas], [temp_likely(end,end) temp_likely(end,:)],[squeeze(std(loglike_all_fun_subfact(end,end,ioff,2:end),[],4)) squeeze(std(loglike_all_fun_subfact(end,:,ioff,2:end),[],4))]);
    
    hold on
    if ioff == noff_all
        title('Sub Fact only')
    end
end
for i = 1:3
    subplot(1,3,i)
    xlim([-22 202])
    xlabel('Orientation')
    ylabel('Log likelihood of 0 deg')
end
suptitle([mouse_list '- Log Likelihood of 0 deg stimulus- all cells- 100X- Sub'])
print(fullfile(LG_base_out, 'loglike_is0_combo_allCells_100X_sub.pdf'),'-dpdf','-fillpage')


figure;
delta_sq_boot = zeros(nDelta,nDelta,noff_all);
delta_sq = zeros(nDelta,nDelta,noff_all);
for ioff = 1:noff_all
    for idelta = 1:nDelta
        for idel = 1:nDelta
            delta_sq_boot(idel,idelta,ioff) = length(find(maxloglike_all(idelta,ioff,2:end) == idel))./1000;
        end
        delta_sq(maxloglike_all(idelta,ioff,1),idelta,ioff) = 1;
    end
    subplot(3,2,(ioff*2)-1)
    imagesc(flipud(delta_sq(:,:,ioff)))
    axis square
    set(gca, 'XTick', 1:nDelta, 'XTickLabels',deltas)
    set(gca, 'YTick', 1:nDelta, 'YTickLabels',fliplr(deltas))
    xlabel('Actual Ori (deg)')
    ylabel('Max Likelihood Ori (deg)')
    title([num2str(chop(off_all(ioff).*(1000/frameRateHz),3)) ' ms ISI'])
    colormap(hot)
    subplot(3,2,(ioff*2))
    imagesc(flipud(delta_sq_boot(:,:,ioff)))
    axis square
    set(gca, 'XTick', 1:nDelta, 'XTickLabels',deltas)
    set(gca, 'YTick', 1:nDelta, 'YTickLabels',fliplr(deltas))
    title([num2str(chop(off_all(ioff).*(1000/frameRateHz),3)) ' ms ISI- bootstrap'])
    xlabel('Actual Ori (deg)')
    ylabel('Max Likelihood Ori (deg)')
    colormap(hot)
    colorbar
end
suptitle([mouse_list '- Max Log Likelihood- All cells- 100X- Sub'])
print(fullfile(LG_base_out, 'maxloglike_byInt_summary_allCells_100X_sub.pdf'),'-dpdf','-fillpage')

%% ROC - sub
h = zeros(nDelta,1);
p = zeros(nDelta,1);
for idelta = 1:nDelta
    [h(idelta,:),p(idelta,:)] = ttest(roc_resp_sub_all(:,1,idelta),roc_resp_sub_all(:,2,idelta));
end
h_roc = zeros(nDelta,nDelta);
p_roc = zeros(nDelta,nDelta);
for i = 1:nDelta
    for idel = 1:nDelta
        ind_cells = intersect(good_ind_theta, find(max_dir_all == idel));
        [h_roc(idel,i), p_roc(idel,i)] = ttest(roc_resp_sub_all(ind_cells,1,i),roc_resp_sub_all(ind_cells,2,i));
    end
end

figure
avg_roc = zeros(2,2,ndiff);
for idiff = 1:ndiff
    del = find(delta_diff ==  diffs(idiff));
    subplot(2,3,idiff)
    scatter(mean(roc_resp_sub_all(good_ind_theta,1,del),3), mean(roc_resp_sub_all(good_ind_theta,2,del),3), 'ob')
    avg_roc(:,1,idiff) = [mean(mean(roc_resp_sub_all(good_ind_theta,1,del),3),1); mean(mean(roc_resp_sub_all(good_ind_theta,2,del),3),1)];
    avg_roc(:,2,idiff) = [std(mean(roc_resp_sub_all(good_ind_theta,1,del),3),[],1)./sqrt(length(good_ind_theta)); std(mean(roc_resp_sub_all(good_ind_theta,2,del),3),[],1)./sqrt(length(good_ind_theta))];
    [htemp ptemp] = ttest(mean(roc_resp_sub_all(good_ind_theta,1,del),3), mean(roc_resp_sub_all(good_ind_theta,2,del),3));
    axis square
    xlim([0 1])
    ylim([0 1])
    refline(1,0)
    vline(0.5)
    hline(0.5)
    xlabel([num2str(chop(off_all(1,:)*(1000/frameRateHz),2)) ' ms ISI'])
    ylabel([num2str(chop(off_all(2,:)*(1000/frameRateHz),2)) ' ms ISI'])
    title(['0 vs ' num2str(diffs(idiff)) '- p = ' num2str(chop(ptemp,2))])
end
subplot(2,3,idiff+1)
errorbar(repmat(diffs, [2 1])', squeeze(avg_roc(:,1,:))', squeeze(avg_roc(:,2,:))', '-o')
ylim([0.3 0.7])
xlim([-10 100])
axis square
ylabel('Average auROC')
xlabel('Diff of Stim from Adaptor (deg)')
suptitle([mouse_list '- ROC 0 (short ISI) vs Diff from Adapter - Sub']) 
print(fullfile(LG_base_out, 'rocv180_allCells_summary_sub.pdf'),'-dpdf','-fillpage')


