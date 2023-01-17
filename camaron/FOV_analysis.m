%% Make subplot of FOVs sorted by depth with title as MOUSE #
% and each figure title as EXP - z_pos. Identify images with most
% non-overlapping cells

clearvars
clc
%% Load experiment info
dataset = 'oriAdapt_V1_cam';
eval(dataset); % run file to load expt.structure

load('Z:\All_Staff\home\camaron\Analysis\2P\good_expt_list.mat')

z_pos = [];


%% Choose list of experiments to view images

% list = i475_good_list;
% list = i472_good_list;
list = i1402_expts;

for iexp = list
    
    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
    CM_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\camaron';
    
    if strcmp(expt(iexp).folder,'lindsey')
        data_base = LG_base;
    elseif strcmp(expt(iexp).folder,'camaron')
        data_base = CM_base;
    end
    
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;

    %% Load image data
    CD = [data_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).runs(1,:)];
    cd(CD);
    imgMatFile = [expt(iexp).runs(1,:) expt(iexp).runs_suffix(1,:) '.mat']; % DONE; Make variable to pull for oriAdapt_V1 that points to imgMatFile of restarted runs (ex: 001_000_001)
    load(imgMatFile);
    z_pos = [z_pos info.config.knobby.pos.z];

end

% sort experiments by depth
z_pos_good_list = [list; z_pos];
z_pos_good_list = flip(sortrows(z_pos_good_list', 2))';


%% Make subplots 

% figure;
[n, n2] = subplotn(length(list));

img_cell = cell(1, length(list));
all_mask_cell = cell(1, length(list));
all_mask_cell_resp = cell(1, length(list));


for i = 1:length(z_pos_good_list)

    iexp = z_pos_good_list(1,i);
    z_pos = z_pos_good_list(2,i);

    run_str = ['runs']; 
    nrun = size(expt(iexp).runs,1);
    for irun = 1:nrun
        run_str = [run_str '-' expt(iexp).runs(irun,:)];
    end

    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    
    %%
    % load redData.mat
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_redData.mat']), 'green_data_avg', 'red_data_avg')
    
    %load mask_cell_red
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']))

    %load responsive cells
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_adaptResp.mat']), 'adapt_resp_ind')
    load(fullfile([data_base '\Analysis\2P'], [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimResp.mat']), 'stim_resp_ind');

    %% Find all responsive cells
    adapt_and_stim_resp_ind = union(adapt_resp_ind{:}, stim_resp_ind);
    resp_cell{i} = adapt_and_stim_resp_ind;

    cells = unique(mask_cell);
    cells = cells(2:end);
    % find values that are NOT in adapt_and_stim_resp_ind, 
    % then set those values to zero
    
    no_resp = cells(~ismember(cells, adapt_and_stim_resp_ind));
    
    mask_cell_resp = mask_cell;
    
    for j = 1:length(no_resp)
        cell_num = no_resp(j);
        mask_cell_resp(mask_cell_resp == cell_num) = 0;
    end
    
    %%
%     imagesc(red_data_avg), colormap gray; caxis([200 1000]);
%     subplot(n,n2,i)
%     sgtitle([mouse ' w/ green__data__avg'])
%     imagesc(green_data_avg)
    img_cell{i} = green_data_avg;
    all_mask_cell{i} = mask_cell;
    all_mask_cell_resp{i} = mask_cell_resp;

%     sgtitle([mouse ' w/ data__dfof__max'] )
%     imagesc(data_dfof_max)
%     title([num2str(iexp) ' @ z= ' num2str(z_pos)])

%     hold on
%     bound = cell2mat(bwboundaries(mask_cell_red(:,:,1))); % 
%     bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
%     plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); hold on;

  
end

%% Create groups of images based on image depth and a given margin of error

depth_list = z_pos_good_list(2,:);

% img_cell

groups = {0} % fill the cell with something to start, remove at the end
exp_ids = z_pos_good_list(1,:)
margin = 25; % for range of z values to create groups

for i = 1:length(z_pos_good_list(2,:)) % for each item in list
    %if any items are between current image and +- margin of this image depth, put into same group
    % if current image is already in a group, skip that image
    if ~ismember(exp_ids(i), groups{end})
        curr_img_depth = z_pos_good_list(2,i);
        tf = (curr_img_depth-margin < z_pos_good_list(2,:)) & (z_pos_good_list(2,:) < curr_img_depth+margin);
        groups{end+1} = exp_ids(tf);
    end
end

groups(1) = []; % remove placeholder value 


%% Once each group is obtained, show groups (with multiple values)
matched_groups = [];
match_inputs = [];

for i = 1:length(groups)

    if length(groups{i}) > 1
        [tf,idx] = ismember(groups{i}, exp_ids);
        grouped_images = img_cell(idx);
        [n n2] = subplotn(length(grouped_images));

        figure();
        
        s = 1; % for subplot indexing 
        for ii = idx
            subplot(n,n2,s)
%             exp_idx = idx(ii)

            imagesc(grouped_images{s})
            title([num2str(exp_ids(ii)) ' @ z= ' num2str(depth_list(ii))])
            s = s+1;
        end
        
        shg % Show current figure
        % option to create prompt so that user can input values to choose
        % images, currently obsolete in this pipline
%         formatSpec = "Create array of %d increasing values for images [%s] (E.g [1 1 2 3]): ";
%         prompt = sprintf(formatSpec, length(groups{i}), num2str(groups{i}));
%         x = input(prompt);
          x = 0; % place holder while testing

        matched_groups = [matched_groups groups(i)];
        match_inputs = [match_inputs {x}];
    end
end

% matched_groups_inputs = [matched_groups ; match_inputs] % row 1 is expt_num, row 2 is numbers for matching expts

cd('Z:\All_Staff\home\camaron\Analysis\2P')
filename = [mouse '_image_matching.mat'];
save(filename, 'mouse', 'z_pos_good_list', 'img_cell', 'all_mask_cell', 'groups', 'margin', 'matched_groups', 'all_mask_cell_resp', 'adapt_and_stim_resp_ind')

close all

%% View groups from saved file - Plot and Save %% Show Masks with significantly responsive cells

% load('Z:\All_Staff\home\camaron\Analysis\2P\i475_image_matching.mat')

exp_ids = z_pos_good_list(1,:);
depth_list = z_pos_good_list(2,:);
matched_group_index = [];

for i = 1:length(groups)

    cd(['Z:\All_Staff\home\camaron\Analysis\2P\' mouse])

    fn = [mouse '_FOV_analysis.pdf'];


    if length(groups{i}) > 1
       matched_group_index = [matched_group_index i];

        [tf,idx] = ismember(groups{i}, exp_ids);
        grouped_images = img_cell(idx);
%         grouped_masks = all_mask_cell(idx); 
        grouped_masks = all_mask_cell_resp(idx); % %% Show Masks with significantly responsive cells
        [n n2] = subplotn(length(grouped_images));



        figure();
        s = 1; % for subplot indexing 
        for ii = idx
            subplot(n,n2,s)

            imagesc(grouped_images{s})

            title(['expt ' num2str(exp_ids(ii)) ' @ z= ' num2str(depth_list(ii))])
            s = s+1;
        end

        exportgraphics(gcf,fn,'Append', true)


        figure(); % Reminder! These masks are for cells that are sig responsive to the adaptor or target?
        s = 1; % for subplot indexing 
        for ii = idx
            subplot(n,n2,s)
%             exp_idx = idx(ii)

%             imagesc(grouped_images{s})
            imagesc(grouped_masks{s})

            title(['expt ' num2str(exp_ids(ii)) ' Resp cells @ z= ' num2str(depth_list(ii))])
            s = s+1;
        end
        
        exportgraphics(gcf,fn,'Append', true)


        figure();
        s = 1; % for subplot indexing 
        for ii = idx
            subplot(n,n2,s)

            shadeimg = imShade(grouped_images{s}, grouped_masks{s});
            imagesc(shadeimg)

            hold on
            curr_mask = grouped_masks{s};
            bound = cell2mat(bwboundaries(curr_mask(:,:,1)));
            plot(bound(:,2),bound(:,1),'.','color','b','MarkerSize',.01); 

            title(['expt ' num2str(exp_ids(ii)) ' @ z= ' num2str(depth_list(ii))])
            s = s+1;
        end

        exportgraphics(gcf,fn,'Append', true)



    end
end

% h = findobj('Type','figure')

%% Next decisions must be made for matching data regarding whether to: keep one/toss others OR match cells and merge trials OR Merge cells across datasets

%Select days with greatest number of non-overlapping cells, filter others
%from groups -> this becomes final expt list 

% Harcoding for now...use algo to acomplish this! % use
% matched_group_index...

% SELECT BEST IMAGE FROM EACH GROUP!
groups_updated = groups;
groups_updated{3} = [120]; % [] value is expt num
groups_updated{4} = [110];
% groups_updated{7} = [72];

% output final expt list

expt_list_final = cell2mat(groups_updated);

cd(['Z:\All_Staff\home\camaron\Analysis\2P\' mouse])
filename = [mouse '_image_matching.mat'];
save(filename, "expt_list_final", '-append')

% dataset = 'oriAdapt_V1_cam';
% eval(dataset); % run file to load expt.structure
% 
% matched_groups = matched_groups_inputs(1,:);
% 
% fprintf('\n')
% 
% 
% for group = 1:length(matched_groups_inputs)
%     list = matched_groups{group};
%     disp(['Group ' num2str(group) ':'])
% 
%     for exp = 1:length(list)
%         iexp = list(exp);
%         disp(['Expt ' num2str(iexp) ': Mouse ' expt(iexp).mouse ' Date: '  expt(iexp).date])
% 
% 
% 
%     end
%         fprintf('\n')
% end

% Check behavior (selectivity) for these days

