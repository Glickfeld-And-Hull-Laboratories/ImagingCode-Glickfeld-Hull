close all
clear all
clc

exp = 3;

ds = 'SST_YM90KDART_trains_ExptList';
eval(ds);

date = expt(exp).date;
abfdate = expt(exp).abfdate;
firstFile = expt(exp).firstFile;
freq_list= expt(exp).freqList;
drug_list = expt(exp).drugList;
omit_list = expt(exp).omit_run_trace;

data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Data\InVitroRecordings\';
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\InVitroRecordings\';
if firstFile<10
    firstFile_str = ['_000' num2str(firstFile)];
else
    firstFile_str = ['_00' num2str(firstFile)];
end
fn_out_use = fullfile(fn_out,date,[abfdate firstFile_str]);

if ~exist(fn_out_use)
    mkdir(fn_out_use)
end

%% all response control and DART
run_ind = find(~isnan(freq_list));
nfile = length(run_ind);
offset = firstFile-run_ind(1);
data = [];

for i = 1:nfile
    run = run_ind(i)+offset;
    if run<10
        [temp si] = abfload(fullfile(data_pn,date,[abfdate '_000' num2str(run) '.abf']));
    else
        [temp si] = abfload(fullfile(data_pn,date,[abfdate '_00' num2str(run) '.abf']));
    end
    if ~isempty(omit_list)
        ind = find(omit_list(:,1) == run_ind(i));
        if ~isempty(ind)
            rem = omit_list(ind,2);
            temp(:,:,rem) = [];
            fprintf(['Remove run ' num2str(run_ind(i)) ' trace ' num2str(rem') '\n']) 
        end
    end
    data(:,i) = squeeze(mean(temp,3));
end

freq_temp = unique(freq_list);
freq = freq_temp(~isnan(freq_temp));
match = zeros(length(freq),2);
for i = 1:length(freq)
    match(i,:) = find(freq_list(run_ind) == freq(i));
end

% resp_start = zeros(nfile,nstim);
% resp_end = zeros(nfile,nstim);
% data_base = zeros(nfile,nstim);
% data_resp = zeros(nfile,nstim);
% 
% figure;
% for i = 1:nfile
%     interval = 1000./freq_list(run_ind(i));
%     subplot(2,3,i)
%     plot(data(:,i))
%     hold on
%     for ii = 1:nstim
%         base_start = firstStim_ind-10+((ii-1)*interval*10);
%         base_end = firstStim_ind+((ii-1)*interval*10);
%         [peak_val peak_ind] = min(data(base_end+20:base_end+60,i));
%         resp_start(i,ii) = base_end+10+peak_ind;
%         resp_end(i,ii) = base_end+30+peak_ind;
%         vline([base_start base_end resp_start(i,ii) resp_end(i,ii)])
%         data_base(i,ii) = mean(data(base_start:base_end,i),1);
%         data_resp(i,ii) = mean(data(resp_start(i,ii):resp_end(i,ii),i),1);
%     end
%     xlim([firstStim_ind 5000])
%     ylim([-500 100])
%     title(num2str(freq_list(run_ind(i))))
% end

resp_start = zeros(nfile,nstim);
resp_end = zeros(nfile,nstim);
data_base = zeros(nfile,nstim);
data_resp = zeros(nfile,nstim);

figure;
for i = 1:nfile
    interval = 1000./freq_list(run_ind(i));
    subplot(2,3,i)
    plot(data(:,i))
    hold on
    for ii = 1:nstim
        base_start = firstStim_ind-10+((ii-1)*interval*10);
        base_end = firstStim_ind+((ii-1)*interval*10);
        if i<4
            [peak_val peak_ind] = min(data(base_end+20:base_end+60,i));
            resp_start(i,ii) = base_end+10+peak_ind;
            resp_end(i,ii) = base_end+30+peak_ind;
        else
            ind = find(match(:,2)==i);
            resp_start(i,ii) = resp_start(match(ind,1),ii);
            resp_end(i,ii) = resp_end(match(ind,1),ii);
        end
        vline([base_start base_end resp_start(i,ii) resp_end(i,ii)])
        data_base(i,ii) = mean(data(base_start:base_end,i),1);
        data_resp(i,ii) = mean(data(resp_start(i,ii):resp_end(i,ii),i),1);
    end
    xlim([firstStim_ind 5000])
    ylim([-500 100])
    title(num2str(freq_list(run_ind(i))))
end

data_sub = data_resp-data_base;
%print(fullfile(fn_out_use,[abfdate firstFile_str '_EPSCTraces.pdf']),'-dpdf','-fillpage')


figure;
for i = 1:size(match,1)
    subplot(size(match,1),1,i)
    for ii = 1:2
        plot(data_sub(match(i,ii),:))
        hold on
    end
    title([num2str(freq_list(match(i,1))) ' Hz'])
    ylim([-600 0])
    xlim([0 11])
    ylabel('EPSC amp (pA)')
end
xlabel('Stimulus #')
subplot(3,1,1)
legend({'Control','YM90K-DART'})
sgtitle('Raw EPSC amplitude')
%print(fullfile(fn_out_use,[abfdate firstFile_str '_rawEPSCAmp.pdf']),'-dpdf','-fillpage')

figure;
for i = 1:size(match,1)
    subplot(size(match,1),1,i)
    for ii = 1:2
        plot(data_sub(match(i,ii),:)./data_sub(match(i,ii),1))
        hold on
    end
    title([num2str(freq_list(match(i,1))) ' Hz'])
    ylim([0 10])
    xlim([0 11])
    ylabel('Norm EPSC')
end
xlabel('Stimulus #')
subplot(3,1,1)
legend({'Control','YM90K-DART'})
sgtitle('Norm to baseline EPSC amplitude')
%print(fullfile(fn_out_use,[abfdate firstFile_str '_normToBaseline.pdf']),'-dpdf','-fillpage')


figure;
for i = 1:size(match,1)
    subplot(size(match,1),1,i)
    for ii = 2
        fprintf([num2str(freq_list(match(i,1))) ' Hz /n'])
        data_sub(match(i,ii),1)./data_sub(match(i,1),1)
        plot(data_sub(match(i,ii),:)./data_sub(match(i,1),:))
        hold on
    end
    title([num2str(freq_list(match(i,1))) ' Hz'])
    ylim([0 1])
    xlim([0 11])
    ylabel('Norm EPSC')
end
xlabel('Stimulus #')
subplot(3,1,1)
sgtitle('Norm to control EPSC amplitude')
%print(fullfile(fn_out_use,[abfdate firstFile_str '_normToControl.pdf']),'-dpdf','-fillpage')

save(fullfile(fn_out_use,[abfdate firstFile_str '_dataSub.mat']),'data', 'data_resp', 'data_base', 'resp_start','resp_end','base_start','base_end','data_sub','match','freq_list')

%% First response control, wash and DART

nfile = length(freq_list);
offset = firstFile-1;
preWashStim = find(drug_list=='w',1,'first');
postWashStim = find(drug_list=='d',1,'first');
data = [];

for i = 1:nfile
    run = i+offset;
    if run<10
        [temp si] = abfload(fullfile(data_pn,date,[abfdate '_000' num2str(run) '.abf']));
    else
        [temp si] = abfload(fullfile(data_pn,date,[abfdate '_00' num2str(run) '.abf']));
    end
    if i == preWashStim
        ind_wash = size(data,2)+1;
    end
    if i == postWashStim
        ind_washend = size(data,2)+1;
    end
    data = [data squeeze(temp)];
end

avg_resp = mean(data,2);

%EPSC across trials
base_start = firstStim_ind-10;
base_end = firstStim_ind;
[peak_val peak_ind] = min(avg_resp(base_end+20:base_end+60,:));
resp_start = base_end+10+peak_ind;
resp_end = base_end+30+peak_ind;
figure; plot(avg_resp)
vline([base_start base_end resp_start resp_end])
data_base = mean(data(base_start:base_end,:),1);
data_resp = mean(data(resp_start:resp_end,:),1);
data_sub = data_resp-data_base;

figure; 
subplot(2,1,1)
plot(data_sub)
hold on
vline([ind_wash ind_washend])
xlabel('Trial')
ylabel('EPSC amplitude')

%Access across trials
[peak_val peak_ind] = min(data(730:740,:),[],1);
base_val = mean(data(720:730,:),1);
access_val = peak_val-base_val;
subplot(2,1,2)
plot(access_val)
hold on
vline([ind_wash ind_washend])
xlabel('Trial')
ylabel('Access current')

print(fullfile(fn_out_use,[abfdate firstFile_str '_respByTrial.pdf']),'-dpdf','-fillpage')
