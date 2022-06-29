clear vars
close all
%%
% Load data

data_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
mouse = 'i475';

data_b = load(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\2P\pooled_data_figs\' mouse '\b_post_processing.mat']);
data_p = load(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\2P\pooled_data_figs\' mouse '\p_post_processing.mat']);
%%
days_b = data_b.data_adapt_dfof_all;
days_p = data_p.data_adapt_dfof_all;

days_b = days_b(~cellfun('isempty',days_b));
days_p = days_p(~cellfun('isempty',days_p));

dates = data_b.dates(~cellfun('isempty',data_b.dates));

resp_cells_days_b = data_b.adapt_resp_ind_cell(~cellfun('isempty',data_b.adapt_resp_ind_cell));
resp_cells_days_p = data_p.adapt_resp_ind_cell(~cellfun('isempty',data_p.adapt_resp_ind_cell));

%% Hubel - Matlab 2022

% cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\2P\pooled_data_figs\i475\i475_resp_amplitude\')


% for iDate = 1:length(dates)
%     curr_day_b = days_b{iDate};
%     curr_day_p = days_p{iDate};
% 
%     % Plot behaving and passive on same figure
% 
%     n_cells = size(curr_day_b,2);
%     [n, n2] = subplotn(n_cells);
% 
%     gca = figure;
% 
%     
%     for i = 1:n_cells
%         cell_b = squeeze(curr_day_b(:,i,:));
%         cell_p = squeeze(curr_day_p(:,i,:));
%         subplot(n,n2,i)
%     
% 
%         shadedErrorBar([], mean(cell_b,2, 'omitnan'), (std(cell_b,0,2, 'omitnan') / sqrt(length(cell_b))),'b');
%         hold on 
%         shadedErrorBar([], mean(cell_p,2, 'omitnan'), (std(cell_p,0,2, 'omitnan') / sqrt(length(cell_p))),'r');
%         vline([20 31 42 53], '--k')
%         vline(64, '--b')  
%         title(['cell ' num2str(i)])
% %         xlabel('Frame')
% %         ylabel('dF/F')
%     end
%     sgtitle(dates{iDate})
% 
% %     print(dates{iDate}, '-dpdf','-fillpage')
%     exportgraphics(gca, 'Z:\All_Staff\home\camaron\Analysis\2P\pooled_data_figs\i475\i475_resp_amplitude\traces_all_cells.pdf', 'Append', true)
% 
% end

%% Nuke - Matlab 2020

cd('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron\Analysis\2P\pooled_data_figs\i475\i475_resp_amplitude\')


for iDate = 1:length(dates)
    curr_day_b = days_b{iDate};
    curr_day_p = days_p{iDate};
    resp_cells_b = cell2mat(resp_cells_days_b{iDate});
    resp_cells_p = cell2mat(resp_cells_days_p{iDate});

    % Plot behaving and passive on same figure

    n_cells = size(curr_day_b,2);
%     [n, n2] = subplotn(n_cells);

    gca = figure;
    sgtitle(dates{iDate})

    
    if n_cells < 20
        [n, n2] = subplotn(n_cells);
    else
        [n, n2] = subplotn(20);
    end
    
%     if n_cells >25
%         nC = 25;
%     else
%         nC = nCells;
%     end
    start = 1;
    x = 0;
    
    for i = 1:n_cells
        if start > 20
            exportgraphics(gca, ['Z:\All_Staff\home\camaron\Analysis\2P\pooled_data_figs\i475\i475_resp_amplitude\' mouse '_traces_all_cells.pdf'], 'Append', true)
            start = 1;
            x = x+1;
            gca = figure;
            sgtitle(dates{iDate})

        end
        

        cell_b = squeeze(curr_day_b(:,i,:));
        cell_p = squeeze(curr_day_p(:,i,:));
        gcs = subplot(n,n2,i-(x.*20));
    
        shadedErrorBar([], mean(cell_b,2, 'omitnan'), (std(cell_b,0,2, 'omitnan') / sqrt(length(cell_b))),'lineProps', 'b');
        hold on 
        shadedErrorBar([], mean(cell_p,2, 'omitnan'), (std(cell_p,0,2, 'omitnan') / sqrt(length(cell_p))),'lineProps', 'r');
        vline([20 31 42 53], '--k')
        vline(64, '--b')  
        
        title(['cell ' num2str(i) '; b:' num2str(ismember(i,resp_cells_b)) '; p:' num2str(ismember(i,resp_cells_p))])
%         xlabel('Frame')
%         ylabel('dF/F')
        set(gcs,'xticklabel',[])
        
        start = start+1;

    end
    


%     print(dates{iDate}, '-dpdf','-fillpage')
    exportgraphics(gca, ['Z:\All_Staff\home\camaron\Analysis\2P\pooled_data_figs\i475\i475_resp_amplitude\' mouse '_traces_all_cells.pdf'], 'Append', true)

end
