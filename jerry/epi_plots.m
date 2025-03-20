clear all; clear global; close all;
clc
ds = 'DART_expt_info'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);

fn_dat_folder = 'G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis';
fn_out = 'G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\plots';

% mkdir(fn_dat_folder);
% mkdir(fn_out);

load('G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\std_cell_PV.mat');
load('G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\std_cell_SST.mat');

%% PV line

% line plot
f1 = figure;
f1.Name = 'PV_std_line';
sgtitle('PV STD Control vs. Post-DART');
xlim([0 3]);
xticks([1 2]);
xticklabels({'Control','DART'});
ylim([0 140]);
hold on

for iMouse = 1:size(std_cell_PV,1)
    curr_mouse = convertCharsToStrings(std_cell_PV{iMouse,2});
    mouse_data = std_cell_PV{iMouse,1};
    if curr_mouse == "i3309"
        plot(mouse_data(1:2),"Marker",'o',"Color","#0047AB","LineWidth",1); % cobalt blue
    elseif curr_mouse == "i3310"
        plot(mouse_data(1:2),"Marker",'o',"Color","#89CFF0","LineWidth",1); % baby blue
    elseif curr_mouse == "i3311"
        plot(mouse_data(1:2),"Marker","^","Color","#D95319","LineWidth",1); % mute orange-ish
    end
end

legend('i3309','i3310','i3311 (PEG)');
hold off

saveas(gcf,fullfile(fn_out,'PV_std_line.pdf'));
%% PV scatter

moi = std_cell_PV(2:4,1);
allMice_pv = cell2mat(moi);
allMice_pv = reshape(allMice_pv,3,3);
allMice_pv = allMice_pv(1:2,:);

f2 = figure;
f2.Name = 'PV_std_scatter';
hold on
scatter(allMice_pv(1,:),allMice_pv(2,:));
sgtitle('PV Scatter Control & DART');
xlim([0 140]);
ylim([0 140]);
ylabel('DART');
xlabel('Control');
plot([0 max([ylim xlim])], [0 max([xlim ylim])], '--r');

hold off
saveas(gcf,fullfile(fn_out,'PV_std_scatter.pdf'));

%% SST plots

% line plot
f3 = figure;
f3.Name = 'SST_std_line';
sgtitle('SST STD Control vs. Post-DART');
xlim([0 3]);
xticks([1 2]);
xticklabels({'Control','DART'});
ylim([0 140]);
hold on

for iMouse = 1:size(std_cell_SST,1)
    curr_mouse = convertCharsToStrings(std_cell_SST{iMouse,2});
    mouse_data = std_cell_SST{iMouse,1};
    if curr_mouse == "i2062"
        plot(mouse_data(1:2),"Marker",'o',"Color","#0047AB","LineWidth",1); % cobalt blue
    elseif curr_mouse == "i2067"
        plot(mouse_data(1:2),"Marker",'o',"Color","#89CFF0","LineWidth",1); % baby blue
    elseif curr_mouse == "i2066"
        plot(mouse_data(1:2),"Marker",'o',"Color","#7393B3","LineWidth",1); % blue grey
    end
end

legend('i2062','i2067','i2066');
hold off

saveas(gcf,fullfile(fn_out,'SST_std_line.pdf'));

% 

%% SST scatter

allMice_sst = cell2mat(std_cell_SST(:,1));

f4 = figure;
f4.Name = 'SST_std_scatter';
hold on
scatter(allMice_sst(:,1),allMice_sst(:,2));
sgtitle('SST Scatter Control & DART');
xlim([0 140]);
ylim([0 140]);
ylabel('DART');
xlabel('Control');
plot([0 max([ylim xlim])], [0 max([xlim ylim])], '--r');

hold off
% 
saveas(gcf,fullfile(fn_out,'SST_std_scatter.pdf'));

%% all scatter

f5 = figure;
f5.Name = 'combine_scatter';
sgtitle('F Scatter');
hold on
scatter(allMice_pv(1,1:2),allMice_pv(2,1:2),"LineWidth",1.5,"MarkerEdgeColor","#0047AB");
scatter(allMice_pv(1,3),allMice_pv(2,3),"LineWidth",1.5,"MarkerEdgeColor","#89CFF0");
scatter(allMice_sst(:,1),allMice_sst(:,2),"LineWidth",1.5,"MarkerEdgeColor","#D95319");
xlim([0 140]);
ylim([0 140]);
plot([0 max([ylim xlim])], [0 max([xlim ylim])], '--',"Color",'black');
legend('PV+DART','PV+PEG','SST+DART');
hold off

saveas(gcf,fullfile(fn_out,'combined_scatter.pdf'));

%% all dots

% first calculate normalized value
norm_pv = allMice_pv(2,:) ./ allMice_pv(1,:);
new_sst = allMice_sst';
norm_sst = new_sst(2,:) ./ new_sst(1,:);

f6 = figure;
f6.Name = 'combined_norm_std';
sgtitle('Normalized STD')
hold on
xticklabels({'PV & DART', 'PV & PEG','SST & DART'});
xticks([1.1 1.7 2.3])
xlim([0.5 2.9])
ylim([0 5])
plot(1.1,norm_pv(1:2),'o','Color',"#0047AB");
plot(1.7,norm_pv(3),'o','Color',"#89CFF0");
plot(2.3,norm_sst,'o','Color',"#D95319");
hold off
ylabel('Post-DART STD/Baseline STD')

% ax = gca;
% ax.FontSize = 8; 

saveas(gcf,fullfile(fn_out,'combined_norm_values.pdf'));