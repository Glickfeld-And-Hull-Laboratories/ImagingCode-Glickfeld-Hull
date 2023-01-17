%% Experiment path
close all
clear all global
clc
dataset = 'oriAdapt_V1_cam';
eval(dataset);
iexp = 48;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
if strcmp(expt(iexp).folder,'lindsey')
    data_base = LG_base;
elseif strcmp(expt(iexp).folder,'camaron')
    data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
end
mouse = expt(iexp).mouse;
date = expt(iexp).date;
pass_str = ['runs-' expt(iexp).pass_run];
%% Load data
load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_adaptResp.mat']))

%% Plot - PASSIVE ONLY

n_cells = size(data_adapt_dfof,2);

gca = figure;
sgtitle(date)

if n_cells < 20
    [n, n2] = subplotn(n_cells);
else
    [n, n2] = subplotn(20);
end

start = 1;
x = 0;

resp_cells = cell2mat(adapt_resp_ind);

for i = 1:n_cells
    if start > 20
%         exportgraphics(gca, ['Z:\All_Staff\home\camaron\Analysis\2P\pooled_data_figs\i475\i475_resp_amplitude\' mouse '_traces_all_cells.pdf'], 'Append', true)
        exportgraphics(gca, fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_traces_all_cells.pdf']), 'Append', true)

        start = 1;
        x = x+1;
        gca = figure;
        sgtitle(date)

    end
    

    cell = squeeze(data_adapt_dfof(:,i,:));
    gcs = subplot(n,n2,i-(x.*20));

%     shadedErrorBar([], mean(cell_b,2, 'omitnan'), (std(cell_b,0,2, 'omitnan') / sqrt(length(cell_b))),'lineProps', 'b');
    hold on 
    shadedErrorBar([], mean(cell,2, 'omitnan'), (std(cell,0,2, 'omitnan') / sqrt(length(cell))),'lineProps', 'r');
    vline([20 31 42 53], '--k')
    vline(64, '--b')  
    
    title(['cell ' num2str(i) '; p:' num2str(ismember(i,resp_cells))])
    xlabel('Frame')
    ylabel('dF/F')
    set(gcs,'xticklabel',[])
    
    start = start+1;

end



% exportgraphics(gca, ['Z:\All_Staff\home\camaron\Analysis\2P\pooled_data_figs\i475\i475_resp_amplitude\' mouse '_traces_all_cells.pdf'], 'Append', true)
exportgraphics(gca, fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' pass_str], [date '_' mouse '_' pass_str '_traces_all_cells.pdf']), 'Append', true)
