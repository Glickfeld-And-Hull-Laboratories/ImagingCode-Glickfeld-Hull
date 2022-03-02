%% Experiment info
clear all; clear global; close all
clc

%Path names
fn_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
fn_analysis = fullfile(fn_base, 'home\Tierney\Analysis\2P'); % CHANGE THIS TO MULTIDAY IMAGING FOLDER
fn_pop = fullfile(fn_base, 'home\Tierney\Analysis\2P\PopulationAnalysis');
% fn_out = fullfile(fn_pop,['population_' time_str,'_',expt(day_id(2)).experiment]);

ds = 'ExperimentData_TD'; % dataset info 
eval(ds)

% day_id_list = 

% Baseline (1), post-MD (2), and recovery (3) days
day_id(2) = 4;
day_id(1) = expt(day_id(2)).matchday_baseline;
day_id(3) = expt(day_id(2)).matchday_recovery;

nd = length(day_id);

%Specific experiment information
mouse = expt(day_id(1)).mouse;

for id = 1:nd
    date = expt(day_id(id)).date;
    runs = expt(day_id(id)).stimruns;
    nrun = length(runs);
    contra = strcmp(expt(day_id(id)).eye_str,'Contra'); % 1 is contra eye open; 0 is ipsi eye open
    
    run_str = catRunName(runs, nrun);
    datemouse = [date '_' mouse];
    datemouserun = [date '_' mouse '_' run_str];
   
    % Load respData and ODI
    fn_respData = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun, '_respData' '.mat']);
    fn_ODI = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun, '_ODI' '.mat']);
    fn_oriResp = fullfile(fn_analysis, datemouse, datemouserun, [datemouserun, '_oriResp' '.mat']);
    temp_respData = load(fn_respData);
    temp_ODI = load(fn_ODI);
    temp_oriResp = load(fn_oriResp);
      
    respData{id} = temp_respData;
    ODI{id} = temp_ODI;
    oriResp{id} = temp_oriResp;
    
    clear temp_respData
    clear temp_ODI
    clear temp_oriResp
end

% if expt(day_id(2)).multiday_time_days>0
%     time_str = [num2str(expt(day_id(2)).multiday_time_days) 'Days'];
% else
%     time_str = 'baseline';
% end


 %% Population averages for ODI
 
 % Get ODIs
 ODI_list = [];
 
 for id = 1:nd
     ODI_list = ODI{id}.ODI;
     mean_ODI(id) = nanmean(ODI_list);
 end
 
swarmchart(day_id(1:3),mean_ODI)
axis square
xlabel('Timepoint')
set(gca,'XTick',3:5,'XTickLabel',{'Pre','Post','Rec'})
ylabel('Mean ODI')
ylim([-1 1])
title('Mean ODI across days')

% MAKE EACH MOUSE A DIFFERENT COLOR

print(fullfile(fn_pop, [mouse '_meanODI.pdf']),'-dpdf','-bestfit')

%% Population average for contra/ipsi inputs

%Get responses from individual eyes
real_ipsi = [];
real_contra = [];
 
 for id = 1:nd
     real_ipsi = ODI{id}.real_ipsi_resp;
     mean_ipsi(id) = mean(real_ipsi);
     real_contra = ODI{id}.real_contra_resp;
     mean_contra(id) = mean(real_contra);
 end

mean_resp = [mean_ipsi; mean_contra]
bar(mean_resp')
axis square
xlabel('Timepoint')
set(gca,'XTick',1:3,'XTickLabel',{'Pre','Post','Rec'})
ylabel('Mean response (df/f)') % CHECK THAT THIS IS ACTUALLY DF/F
legend('Ipsi','Contra')
title(['Ipsi v. contra response'])

print(fullfile(fn_pop, [mouse '_EyeResp.pdf']),'-dpdf','-bestfit')

save(fullfile(fn_out, ['PopulationData.mat']), 'mean_ODI', 'mean_ipsi', 'mean_contra')
%% Population average for tuning.

Eyes = unique(contra);
nEye = length(Eyes);

%Get tuning widths (k from von mises) from individual eyes
figure
titles = {'Ipsi', 'Contra'} % CHECK THAT THIS IS THE RIGHT ORDER. REMINDER HOW TO NOT HARD CODE???
for iEye = 1:nEye
    subplot(1,2,iEye)
    for id = 1:nd
        cdfplot(oriResp{id}.k1_ori(iEye,:));
        hold on
    end
    axis square
    title(titles{iEye})
end
hold off
legend('Pre','Post','Rec','Location','Southeast');

print(fullfile(fn_pop, [mouse '_TuningCDF.pdf']),'-dpdf','-bestfit')

% Ipsi eye -- MD v. baseline
subplot(2,2,1)
scatter(oriResp{1}.k1_ori(1,:),oriResp{2}.k1_ori(1,:));

% Ipsi eye -- MD v. recovery
subplot(2,2,2)
scatter(oriResp{1}.k1_ori(1,:),oriResp{3}.k1_ori(1,:));

% Contra eye -- MD v. baseline
subplot(2,2,3)
scatter(oriResp{1}.k1_ori(1,:),oriResp{2}.k1_ori(2,:));

% Contra eye -- MD v. recovery
subplot(2,2,4)
scatter(oriResp{1}.k1_ori(1,:),oriResp{3}.k1_ori(2,:));

print(fullfile(fn_pop, [mouse '_TuningScatter.pdf']),'-dpdf','-bestfit')