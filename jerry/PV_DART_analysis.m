clear all; clear global; close all;
clc

ds = 'DART_expt_info'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);

% day_id = [24 29 32 35 38 40 42 44]; % concat days, enter one post-DART day id for each mouse

fnroot = fullfile(rc.achAnalysis,'PV_YM90K','summary_analyses');
fn_epi = fullfile(fnroot,'epileptiform');
cd(fn_epi);

nOn = 30;
nOff = 60;
nTot = nOn+nOff;
ntrials = 960;

%%
load('prepped_np_TCs.mat');
%% PV_DART off_np_tc epileptiform quantification

% off_np_tc: each row is a mouse. Each mouse cell contains a n x 2 cell array
% where n is the number of imaging sessions for the mouse, the 1st column
% is an averaged-across-cells timecourse made from a concatenation of the last 30 off-frames of
% every trial, and the 2nd column is the mean-subtracted timecourse of the
% 1st column. 

cd(fn_epi);
load('off_np_tc.mat');
nMice = size(off_np_tc,1);
perc_mat = NaN(nMice,3);
std_mat = NaN(nMice,3);
t_since_drug_mat = strings(nMice,3);

for mouse = 1:nMice %iterate through mouse
    mouse_data = off_np_tc{mouse,1};
    mouse_id = off_np_tc{mouse,2};
    this_sessions = query_expt(str2double(mouse_id(2:end)));

    all_raw_f = cat(1,mouse_data{:,1});
    ylim_max = 20 + max(all_raw_f);
    ylim_min = min(all_raw_f)-20;

    figure;
    sgtitle(mouse_id);
    for sess_idx = 1:size(mouse_data,1)
        curr_sesh = mouse_data{sess_idx,1}; % data for the current imaging session of the current mouse
        this_session_num = this_sessions(sess_idx);
        t_since_drug_hrs = expt(this_session_num).multiday_timesincedrug_hours;
        t_since_drug_mat(mouse,sess_idx) = t_since_drug_hrs;
        nplots = length(this_sessions);

        if t_since_drug_hrs == '0'
            this_tp_title = 'baseline';
        else
            this_tp_title = [t_since_drug_hrs ' hrs post-DART'];
        end

        mean_curr_sesh = mean(curr_sesh);
        std_curr_sesh = std(curr_sesh);
        frac_past = sum(curr_sesh > (mean_curr_sesh+3*std_curr_sesh))/length(curr_sesh);
        perc_mat(mouse, sess_idx) = frac_past;
        std_mat(mouse, sess_idx) = std_curr_sesh;

        subplot(nplots,1,sess_idx);
        plot(curr_sesh);
        ylim([ylim_min ylim_max]);
        title([this_tp_title ' std= ' num2str(std_curr_sesh) ' frac=' num2str(frac_past)]);
    end
    % saveas(gcf,fullfile(fn_epi,'plots',[mouse_id '.pdf']));
end

save(fullfile(fn_epi,'epi_data.mat'),"std_mat","perc_mat","t_since_drug_mat");
clear all_raw_f ylim_max ylim_min std_mat perc_mat curr_sesh off_np_tc
%% plot epilepsy data as modulation index

cd(fn_epi);
load('epi_data.mat');
load('DART_expt_info.mat');
% for dart_expt_info, 2nd column is dart dosage, 3rd column is virus titer
% DART:
% 1 = YM90K-PEG 10:1 alx647-COOH
% 2 = 1:10:1 YM90K-DART:blank-DART:alx-DART
% 3 = 4:6:1 YM90K-DART:blank-DART:alx-DART
% 4 = 5:5:1 YM90K-DART:blank-DART:alx-DART
% 5 = YM90K-DART 10:1 alx647-DART
% virus:
% 1 = GCaMP6s 1:1 MB HTP
% 2 = GCaMP8f 1:1 NES HTP
% 3 = 2:1:1 GCaMP8f:NES HTP:aCSF

if exist('off_np_tc','var') == 1
    mouse_list = convertCharsToStrings(off_np_tc(:,2)); 
    clear off_np_tc
else
    load('off_np_tc.mat');
    mouse_list = convertCharsToStrings(off_np_tc(:,2));
    clear off_np_tc
end

nMice = size(mouse_list,1);
TPs_mat = t_since_drug_mat(:,2:3);
std_index_mat = nan(nMice,3);
frac_index_mat = nan(nMice,3);

% Fig 1 frac norm
figure;
set(gcf, 'Name', 'EpiPlot_frac_norm');
title('Modulation of Epileptic Activity, frac, normalized');
xlim([0 26]);
ylim([0 3.7]);
xlabel('Time Since DART Infusion (hrs)');
ylabel('Fraction Modulation Index (AU)');
hold on

for mouse = 1:nMice
    if sum(~isnan(std_mat(mouse,:))) == 2
        nTP = 2;
    elseif sum(~isnan(std_mat(mouse,:))) == 3
        nTP = 3;
    else
        error('Incorrect number of timepoints, check data matrix dimensions')
    end
    % calculate normalized fraction value
    for tp = 1:nTP
        std_index_mat(mouse,tp) = std_mat(mouse,tp)/std_mat(mouse,1);
        frac_index_mat(mouse,tp) = perc_mat(mouse,tp)/perc_mat(mouse,1);
    end

    % color code by dart dosage 
    if expt_info{mouse,2} == 1
        this_col = "#D95319"; % mute orange-ish
    elseif expt_info{mouse,2} == 2
        this_col = "#40E0D0"; % turquoise
    elseif expt_info{mouse,2} == 3
        this_col = "#89CFF0"; % baby blue
    elseif expt_info{mouse,2} == 4
        this_col = "#7393B3"; % blue gray
    elseif expt_info{mouse,2} == 5
        this_col = "#0047AB"; % cobalt blue
    else
        error('Unrecognized DART dosage value. Check metadata values.')
    end

    % line style code by virus titer/type
    if expt_info{mouse,3} == 1
        this_marker = "o"; 
    elseif expt_info{mouse,3} == 2
        this_marker = "^"; 
    elseif expt_info{mouse,3} == 3
        this_marker = "square"; 
    else
        error('Unrecognized viral titer. Check metadata values.');
    end
    plot(str2double(t_since_drug_mat(mouse,:)),frac_index_mat(mouse,:),"Color",this_col,"Marker",this_marker,"LineWidth",1.5);
end

legend
hold off

% fig 2 frac un-norm

figure;
set(gcf, 'Name', 'EpiPlot_frac_raw');
title('Modulation of Epileptic Activity, frac, raw value');
xlim([0 26]);
ylim([0 0.027]);
xlabel('Time Since DART Infusion (hrs)');
ylabel('Fraction of Frames');
hold on

for mouse = 1:nMice
    if sum(~isnan(std_mat(mouse,:))) == 2
        nTP = 2;
    elseif sum(~isnan(std_mat(mouse,:))) == 3
        nTP = 3;
    else
        error('Incorrect number of timepoints, check data matrix dimensions')
    end

    % color code by dart dosage 
    if expt_info{mouse,2} == 1
        this_col = "#D95319"; % mute orange-ish
    elseif expt_info{mouse,2} == 2
        this_col = "#40E0D0"; % turquoise
    elseif expt_info{mouse,2} == 3
        this_col = "#89CFF0"; % baby blue
    elseif expt_info{mouse,2} == 4
        this_col = "#7393B3"; % blue gray
    elseif expt_info{mouse,2} == 5
        this_col = "#0047AB"; % cobalt blue
    else
        error('Unrecognized DART dosage value. Check metadata values.')
    end

    % line style code by virus titer/type
    if expt_info{mouse,3} == 1
        this_marker = "o"; 
    elseif expt_info{mouse,3} == 2
        this_marker = "^"; 
    elseif expt_info{mouse,3} == 3
        this_marker = "square"; 
    else
        error('Unrecognized viral titer. Check metadata values.');
    end
    plot(str2double(t_since_drug_mat(mouse,:)),perc_mat(mouse,:),"Color",this_col,"Marker",this_marker,"LineWidth",1.5);
end
legend
hold off

% Fig 3 std norm
figure;
set(gcf, 'Name', 'EpiPlot_std_norm');
title('Modulation of Epileptic Activity, std, normalized');
xlim([0 26]);
% ylim([0 3.7]);
xlabel('Time Since DART Infusion (hrs)');
ylabel('Std Modulation Index (AU)');
hold on

for mouse = 1:nMice
    if sum(~isnan(std_mat(mouse,:))) == 2
        nTP = 2;
    elseif sum(~isnan(std_mat(mouse,:))) == 3
        nTP = 3;
    else
        error('Incorrect number of timepoints, check data matrix dimensions')
    end
    % calculate normalized fraction value
    for tp = 1:nTP
        std_index_mat(mouse,tp) = std_mat(mouse,tp)/std_mat(mouse,1);
        frac_index_mat(mouse,tp) = perc_mat(mouse,tp)/perc_mat(mouse,1);
    end

    % color code by dart dosage 
    if expt_info{mouse,2} == 1
        this_col = "#D95319"; % mute orange-ish
    elseif expt_info{mouse,2} == 2
        this_col = "#40E0D0"; % turquoise
    elseif expt_info{mouse,2} == 3
        this_col = "#89CFF0"; % baby blue
    elseif expt_info{mouse,2} == 4
        this_col = "#7393B3"; % blue gray
    elseif expt_info{mouse,2} == 5
        this_col = "#0047AB"; % cobalt blue
    else
        error('Unrecognized DART dosage value. Check metadata values.')
    end

    % line style code by virus titer/type
    if expt_info{mouse,3} == 1
        this_marker = "o"; 
    elseif expt_info{mouse,3} == 2
        this_marker = "^"; 
    elseif expt_info{mouse,3} == 3
        this_marker = "square"; 
    else
        error('Unrecognized viral titer. Check metadata values.');
    end
    plot(str2double(t_since_drug_mat(mouse,:)),std_index_mat(mouse,:),"Color",this_col,"Marker",this_marker,"LineWidth",1.5);
end

legend
hold off

% fig 4 std un-norm

figure;
set(gcf, 'Name', 'EpiPlot_std_raw');
title('Modulation of Epileptic Activity, std, raw value');
xlim([0 26]);
% ylim([0 0.027]);
xlabel('Time Since DART Infusion (hrs)');
ylabel('Standard Deviation of F');
hold on

for mouse = 1:nMice
    if sum(~isnan(std_mat(mouse,:))) == 2
        nTP = 2;
    elseif sum(~isnan(std_mat(mouse,:))) == 3
        nTP = 3;
    else
        error('Incorrect number of timepoints, check data matrix dimensions')
    end

    % color code by dart dosage 
    if expt_info{mouse,2} == 1
        this_col = "#D95319"; % mute orange-ish
    elseif expt_info{mouse,2} == 2
        this_col = "#40E0D0"; % turquoise
    elseif expt_info{mouse,2} == 3
        this_col = "#89CFF0"; % baby blue
    elseif expt_info{mouse,2} == 4
        this_col = "#7393B3"; % blue gray
    elseif expt_info{mouse,2} == 5
        this_col = "#0047AB"; % cobalt blue
    else
        error('Unrecognized DART dosage value. Check metadata values.')
    end

    % line style code by virus titer/type
    if expt_info{mouse,3} == 1
        this_marker = "o"; 
    elseif expt_info{mouse,3} == 2
        this_marker = "^"; 
    elseif expt_info{mouse,3} == 3
        this_marker = "square"; 
    else
        error('Unrecognized viral titer. Check metadata values.');
    end
    plot(str2double(t_since_drug_mat(mouse,:)),std_mat(mouse,:),"Color",this_col,"Marker",this_marker,"LineWidth",1.5);
end
legend
hold off

% figHandles = findall(0,'Type','figure');
% this_figname = get(figHandles,'Name');
% for i = 1:numel(figHandles)
%     t = this_figname{i,1};
%     saveas(figHandles(i),fullfile(fn_epi,'plots',[t '.pdf']));
% end
%
% or use "save_all_open_figs(directory_to_be_saved_to)"


%% power frequency curve

% extract just i3309 and i3310, the two high dose YM90K mice
% load('prepped_np_TCs.mat');
% first column stores raw F TC (but only the last 30 off frames of every trial), 
% second stores mean-subtracted (same structure as the first), third stores nan-padded (full trial
% length), fourth is full np tc (full trial length)

mus = prepped_np_TCs(2:4,:);
nMice = size(mus,1);
close all
mus_rawPwr = cell(size(mus,1),2);
mus_rawPwr(:,2) = mus(:,2);
% norm_dataset = cell(nMice,3);
% norm_pks = zeros(nMice,3);
for iMouse = 1:nMice
    mouse = mus{iMouse,1};
    intact_tc = mouse(:,1);
    meanSub_tc = mouse(:,2);
    nanpad_tc = mouse(:,3);
    full_tc = mouse(:,4);
    meanSub_full_tc = mouse(:,5);
    nSesh = size(mouse,1);
    curr_mus = mus{iMouse,2};
    sesh_fits = cell(nSesh,1);
    sesh_db = cell(nSesh,1);
    sesh_f = cell(size(mouse,1),1);
    for sesh = 1:nSesh
        f1 = figure;
        [p,f] = pspectrum(meanSub_tc{sesh},15); 
        dbp = pow2db(p);
        plot(f,p);
        % ylim([-15 30]);
        %xlim([0 100]);
        sgtitle([curr_mus ' Session ' num2str(sesh) ' Power Spectrum'])
        xlabel('frequency (Hz)')
        ylabel('Power')
        f1.Name = [curr_mus '_Session' num2str(sesh) '_ps_noLog'];
        sesh_f{sesh} = p;
        % full sesh PS
        % f1 = figure;
        % [p,f] = pspectrum(meanSub_full_tc{sesh},15); 
        % dbp2 = pow2db(p);
        % plot(f,dbp2);
        % dbp_fit = smooth(dbp2,0.03,'lowess'); 
        % hold on
        % plot(f,dbp_fit,"LineWidth",2);
        % hold off
        % ylim([-15 30]);
        % %xlim([0 100]);
        % sgtitle([curr_mus ' Session ' num2str(sesh) ' Full Session'])
        % xlabel('frequency (Hz)')
        % ylabel('Power (dB)')
        % f1.Name = [curr_mus '_Session' num2str(sesh) '_ps_full_fit'];
        % sesh_fits{sesh} = dbp_fit;
        % sesh_db{sesh} = f;
        
        % f2 = figure;
        % spectrogram(meanSub_full_tc{sesh},900,[],[],15,'yaxis',MinThreshold=-20);
        % clim([-20 55]);
        % sgtitle([curr_mus ' Session ' num2str(sesh) ' Full Session Spectrogram']);
        % f2.Name = [curr_mus '_Session' num2str(sesh) '_TimeFreq'];

        % this_norm_p = normalize(dbp,'range',[0 1]);
        % norm_dataset{iMouse,sesh} = this_norm_p;
        % [norm_pks(iMouse,sesh),idx_max] = max(this_norm_p(2501:end));
        % linex = f(2500+idx_max);
        % 
        % f2 = figure;
        % plot(f,this_norm_p);
        % hold on
        % xline(linex);
        % hold off
        % f2.Name = [curr_mus '_Session' num2str(sesh) '_normalized'];
        % xlabel('frequency (Hz)')
        % ylabel('z-score')
        % sgtitle([curr_mus ' Session ' num2str(sesh) ' Normalized PS']);
        % f2 = figure;
        % plot(meanSub_tc{sesh});
        % ylim([-200 1100]);
        % sgtitle([curr_mus ' Session ' num2str(sesh) ' Mean-Subtracted Neuropil TC']);
        % xlabel('nFrame')
        % ylabel('F')
        % f2.Name = [curr_mus '_Session' num2str(sesh) '_MSNPTC'];
    end
    mus_rawPwr{iMouse,1} = sesh_f;
    % for sesh = 1:nSesh
    %     if sesh == 1
    %         continue
    %     end
    %     this_curve = sesh_fits{sesh} - sesh_fits{1};
    %     this_f = sesh_db{sesh};
    %     f3 = figure;
    %     plot(this_f,this_curve);
    %     ylim([-8 16]);
    %     %xlim([0 100]);
    %     sgtitle([curr_mus ' Session ' num2str(sesh) ' Over Baseline, Full Session'])
    %     xlabel('frequency (Hz)')
    %     ylabel('Power (dB)')
    %     f3.Name = [curr_mus '_Session' num2str(sesh) '_ps_full_fit_ovb'];
    % end
    % pause
    % disp('press any key to continue');
end



% figHandles = findall(0,'Type','figure'); 
% this_fig = get(figHandles,'Name');
% [a,b] = subplotn(size(this_fig,1));
% figure;
% 
% tile_plot = tiledlayout(a,b);
% this_fig = flip(this_fig);
% for i = 1:a
%     for j = 1:b
%        subplot() 
%     end
% end



% save(fullfile(fn_epi,'rawPwrTCs'),'mus_rawPwr','-v7.3');
%save_all_open_figs('G:\home\ACh\Analysis\2p_analysis\PV_YM90K\summary_analyses\epileptiform\plots\power');

%% integral
if exist('mus_rawPwr','var') == 0
    load(fullfile(fn_epi,"rawPwrTCs.mat"))
end

nMice = size(mus_rawPwr,1);

spacing = 7.5/4096;
idx_freq_range = spacing:spacing:7.5;
intValues1 = nan(nMice,3);
intValues2 = nan(nMice,3);
intValues3 = nan(nMice,3);
intValues4 = nan(nMice,3);

for mouse = 1:nMice
    thisMousePs = mus_rawPwr{mouse,1};
    nSesh = length(thisMousePs);
    for sesh = 1:nSesh
        thisCurve = thisMousePs{sesh};
        % integral of curve where freq >= 0.5 Hz
        idx4int = idx_freq_range >= 0.5;
        intValues1(mouse,sesh) = trapz(thisCurve(idx4int));
        % integral of curve across whole freq
        intValues2(mouse,sesh) = trapz(thisCurve);
        % integral of curve where freq <= 0.5 Hz
        idx4int = idx_freq_range <= 0.5;
        intValues3(mouse,sesh) = trapz(thisCurve(idx4int));
        % integral of curve where freq <= 1 Hz
        idx4int = idx_freq_range <= 1;
        intValues4(mouse,sesh) = trapz(thisCurve(idx4int));
    end
end

% T1 = array2table(intValues1,"VariableNames",{'session 1','session 2','session 3'},"RowNames",{'i3309','i3310','i3311'});
% fig1 = uifigure;
% uitable(fig1,"Data",T1,"Position",[0 0 1 1],'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
% lbl1 = uilabel(fig1, 'Position',[100 100 125 250]);
% %lbl1.Text = '$$\int_{0.5}^{7.5} f(freq)\;d freq$$';
% lbl1.Text = "freq $\ge$ 0.5 Hz";
% lbl1.Interpreter = "latex";
% fig1.Name = ('Integral of Power Spectrum, freq >= 0.5 Hz');
% 
% T2 = array2table(intValues2,"VariableNames",{'session 1','session 2','session 3'},"RowNames",{'i3309','i3310','i3311'});
% fig2 = uifigure;
% uitable(fig2,"Data",T2,"Position",[0 0 1 1],'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
% lbl2 = uilabel(fig2, 'Position',[100 100 125 250]);
% %lbl1.Text = '$$\int_{0.5}^{7.5} f(freq)\;d freq$$';
% lbl2.Text = "All freq";
% fig2.Name = ('Integral of Power Spectrum, All freq');
% 
% T3 = array2table(intValues3,"VariableNames",{'session 1','session 2','session 3'},"RowNames",{'i3309','i3310','i3311'});
% fig3 = uifigure;
% uitable(fig3,"Data",T3,"Position",[0 0 1 1],'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
% lbl3 = uilabel(fig3, 'Position',[100 100 125 250]);
% %lbl1.Text = '$$\int_{0.5}^{7.5} f(freq)\;d freq$$';
% lbl3.Text = "freq $\le$ 0.5 Hz";
% lbl3.Interpreter = "latex";
% fig3.Name = ('Integral of Power Spectrum, freq <= 0.5 Hz');
% 
% T4 = array2table(intValues3,"VariableNames",{'session 1','session 2','session 3'},"RowNames",{'i3309','i3310','i3311'});
% fig4 = uifigure;
% uitable(fig4,"Data",T4,"Position",[0 0 1 1],'Units', 'Normalized', 'Position',[0, 0, 1, 1]);
% lbl4 = uilabel(fig4, 'Position',[100 100 125 250]);
% %lbl1.Text = '$$\int_{0.5}^{7.5} f(freq)\;d freq$$';
% lbl4.Text = 'freq $\le$ 1 Hz';
% lbl4.Interpreter = "latex";
% fig4.Name = ('Integral of Power Spectrum, freq <= 1 Hz');

%%
fig1 = figure;
fig1.Name = 'IntPS_fGE0.5Hz_norm';
hold on
% for mouse = 1:nMice
%     plot(log10(intValues1(mouse,:)),"Marker",'o',"LineWidth",1.5);
% end
for mouse = 1:nMice
    plot(intValues1(mouse,:)/intValues1(mouse,1),"Marker",'o',"LineWidth",1.5);
end
hold off
legend('i3309','i3310','i3311');
sgtitle('Normalized Integral of Power Spectrum, freq $\ge$ 0.5 Hz','Interpreter', 'LaTeX')
ylabel('Normalized Integral Value (AU)')
xlabel('Session #')
xticks([1 2 3])
xticklabels({'1','2','3'})

fig2 = figure;
fig2.Name = 'IntPS_FR_norm';
hold on
% for mouse = 1:nMice
%     plot(log10(intValues2(mouse,:)),"Marker",'o',"LineWidth",1.5);
% end
for mouse = 1:nMice
    plot(intValues2(mouse,:)/intValues2(mouse,1),"Marker",'o',"LineWidth",1.5);
end
hold off
legend('i3309','i3310','i3311');
sgtitle('Normalized Integral of Power Spectrum, Full Range','Interpreter', 'LaTeX')
ylabel('Normalized Integral Value (AU)')
xlabel('Session #')
xticks([1 2 3])
xticklabels({'1','2','3'})

fig3 = figure;
fig3.Name = 'IntPS_fLE0.5Hz_norm';
hold on
% for mouse = 1:nMice
%     plot(log10(intValues3(mouse,:)),"Marker",'o',"LineWidth",1.5);
% end
for mouse = 1:nMice
    plot(intValues3(mouse,:)/intValues3(mouse,1),"Marker",'o',"LineWidth",1.5);
end
hold off
legend('i3309','i3310','i3311');
sgtitle('Normalized Integral of Power Spectrum, freq $\le$ 0.5 Hz','Interpreter', 'LaTeX')
ylabel('Normalized Integral Value (AU)')
xlabel('Session #')
xticks([1 2 3])
xticklabels({'1','2','3'})

fig4 = figure;
fig4.Name = 'IntPS_fLE1Hz_norm';
hold on
% for mouse = 1:nMice
%     plot(log10(intValues4(mouse,:)),"Marker",'o',"LineWidth",1.5);
% end
for mouse = 1:nMice
    plot(intValues4(mouse,:)/intValues4(mouse,1),"Marker",'o',"LineWidth",1.5);
end
hold off
legend('i3309','i3310','i3311');
sgtitle('Normalized Integral of Power Spectrum, freq $\le$ 1 Hz','Interpreter', 'LaTeX')
ylabel('Normalized Integral Value (AU)')
xlabel('Session #')
xticks([1 2 3])
xticklabels({'1','2','3'})

save_all_open_figs('G:\home\ACh\Analysis\2p_analysis\PV_YM90K\summary_analyses\epileptiform\plots');
%% manually calculate ps
% Fs = 15; % freq
% T = 1/Fs; % period length
% L = length(meanSub_tc{sesh});  % length of sample       
% t = (0:L-1)*T;  % Time vector
% y = fft(meanSub_tc{sesh});
% figure;
% plot(Fs/L*(0:L-1),abs(y));
% 
% P2 = abs(y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% f = Fs/L*(0:(L/2));
% figure
% plot(f,pow2db(P1.^2),'LineWidth',0.1);
% ylim([-15 25]);

% x = spectrogram(meanSub_tc{sesh});
% f3 = figure;
% plot(x(:,1));
% 
% N = 128;
% x = [1 1/sqrt(2)].*exp(1j*pi./[4;2]*(0:N-1)).';
% 
% [p,f] = pspectrum(x);

% plot(f/pi,p)
% hold on
% stem(0:2/N:2-1/N,abs(fft(x)/N).^2)
% hold off
% axis([0.15 0.6 0 1.1])
% legend("Channel "+[1;2]+", "+["pspectrum" "fft"])
% grid
% 
% Fs = 1000;
% t = (0:1/Fs:0.296)';
% x = cos(2*pi*t*200)+0.1*randn(size(t));
% xTable = timetable(seconds(t),x);
% 
% [pxx,f] = pspectrum(xTable);
% 
% plot(f,pow2db(pxx))
% grid on
% xlabel('Frequency (Hz)')
% ylabel('Power Spectrum (dB)')
% title('Default Frequency Resolution')


%% pwelch

t = 1/15:1/15:86400/15; % Time vector
x = meanSub_tc{sesh}; % Signal 

% Compute the power spectral density using pwelch
[pxx, f] = pwelch(x, [], [], [], 1/(t(2)-t(1)));

% Plot the power spectral density
figure;
plot(f, 10*log10(pxx)); % Plot in dB scal
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density using pwelch');



