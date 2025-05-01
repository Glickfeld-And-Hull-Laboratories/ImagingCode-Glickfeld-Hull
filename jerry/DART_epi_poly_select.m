clear all; clear global; close all;
clc

dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories


fn_analysis_root = 'G:\home\ACh\Analysis\2p_analysis';
fn_data_root = 'G:\home\ACh\Data\2p_data';
fn_out_root = 'G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\polymask';

session_id_TH = [28 31 34 69 71]; % enter post-DART session IDs
session_id_CC = [169 177 183]; % enter post-DART session IDs

%% load data and reg info TH


ds = 'DART_expt_info';
eval(ds);

polyTCcell = cell(length(session_id_TH),2);

for id = 1:length(session_id_TH)
    polyTCcell{id,2} =  expt(session_id_TH(id)).mouse;
    TCmat = cell(2,1);
    pre_day = expt(session_id_TH(id)).multiday_matchdays;
    post_day = session_id_TH(id);
    
    days = [pre_day post_day];
    
    for iday = 1:length(days)
        % load outs
        day = days(iday);
        if iday == 1
            pre_or_post = 'pre';
        elseif iday == 2
            pre_or_post = 'post';
        end

        fn_ana_full = fullfile(fn_analysis_root,expt(day).exptType,expt(day).mouse,expt(day).date,expt(day).contrastxori_runs);
        load(fullfile(cell2mat(fn_ana_full),'regOuts&Img.mat'));
        % find original data
        fn_data_full = fullfile(fn_data_root,expt(day).mouse,expt(day).date,expt(day).contrastxori_runs);
        cd(cell2mat(fn_data_full));
        sesh2load = [cell2mat(expt(day).contrastxori_runs) '_000_000'];
        data_g = sbxread(sesh2load,0,86400);
        data_g = squeeze(data_g(1,:,:,:));
        clear global;
        [~,data_g_reg] = stackRegister_MA(data_g,[],[],outs);
        clear data_g
        % outs_gpu = gpuArray(outs(1:2400,:));
        % data_gpu = gpuArray(data_g);
        % [~,data_gpu_reg] = stackRegister_TH(data_gpu,[],[],outs_gpu);

        disp(cell2mat(['data is from ' expt(day).mouse ' ' expt(day).date ' ' expt(day).contrastxori_runs]));
        % data_g_reg = gather(data_gpu_reg);
        % clear data_gpu_reg outs_gpu data_gpu
        % gpuDevice([]);
        
        avg_img = mean(data_g_reg,3);
        f1=figure;
        imagesc(avg_img);
        f1.WindowState = "maximized";
        title(['Select ROI ' expt(day).mouse ' ' pre_or_post]);
        mask_poly = roipoly();
        close();

        roiPixels = avg_img;
        roiPixels(~mask_poly) = 0;
        f2 = figure;
        f2.Name = ['PolyROI_' expt(day).mouse '_' pre_or_post];
        imagesc(roiPixels);
        title(['Polygon ROI ' expt(day).mouse ' ' pre_or_post]);
        
        TC = stackGetTimeCourses(data_g_reg,mask_poly);
        TCmat{iday,1} = TC;
    end
    polyTCcell{id,1} = TCmat;
end


TCs_PV = polyTCcell;

save(fullfile(fn_out_root,'TCs_PV'),'TCs_PV','-v7.3');


%% 
load(fullfile(fn_out_root,'TCs_PV'));
for iMouse = 1:size(TCs_PV,1)
    this_mouse_data = TCs_PV{iMouse,1};
    mouse_id = TCs_PV{iMouse,2};
    for row = 1:length(this_mouse_data)
        this_mouse_data{row,1} = this_mouse_data{row,1} - mean(this_mouse_data{row,1});
    end

    both_day_data = cell2mat(this_mouse_data);
    ymax = max(both_day_data);
    ymin = min(both_day_data);

    for iday = 1:length(this_mouse_data)
        this_TC = this_mouse_data{iday};
        if iday == 1
            pre_or_post = 'pre';
        elseif iday == 2
            pre_or_post = 'post';
        end

        f1 = figure;
        plot(this_TC)
        f1.Name = ['PolyTC_' mouse_id '_' pre_or_post];
        sgtitle(['Polygon TC ' mouse_id ' ' pre_or_post]);
        ylim([ymin ymax])
        ylabel('F');
        xlabel('Frame');
    end
end

save_all_open_figs(fullfile(fn_out_root,'plots'));

close all
%% load data and reg info CC


ds = 'DART_V1_contrast_ori_Celine';
eval(ds);

polyTCcell = cell(length(session_id_CC),2);

for id = 1:length(session_id_CC)
    polyTCcell{id,2} =  expt(session_id_CC(id)).mouse;
    TCmat = cell(2,1);
    pre_day = expt(session_id_CC(id)).multiday_matchdays;
    post_day = session_id_CC(id);
    
    days = [pre_day post_day];
    
    for iday = 1:length(days)
        % load outs
        day = days(iday);
        if iday == 1
            pre_or_post = 'pre';
        elseif iday == 2
            pre_or_post = 'post';
        end

        fn_ana_full = fullfile(fn_analysis_root,expt(day).exptType,expt(day).mouse,expt(day).date,expt(day).contrastxori_runs);
        load(fullfile(cell2mat(fn_ana_full),'regOuts&Img.mat'));
        % find original data
        fn_data_full = fullfile(fn_data_root,expt(day).mouse,expt(day).date,expt(day).contrastxori_runs);
        cd(cell2mat(fn_data_full));
        sesh2load = [cell2mat(expt(day).contrastxori_runs) '_000_000'];
        data_g = sbxread(sesh2load,0,43200);
        data_g = squeeze(data_g(1,:,:,:));
        clear global;
        [~,data_g_reg] = stackRegister_MA(data_g,[],[],outs);
        clear data_g
        % outs_gpu = gpuArray(outs(1:2400,:));
        % data_gpu = gpuArray(data_g);
        % [~,data_gpu_reg] = stackRegister_CC(data_gpu,[],[],outs_gpu);

        disp(cell2mat(['data is from ' expt(day).mouse ' ' expt(day).date ' ' expt(day).contrastxori_runs]));
        % data_g_reg = gather(data_gpu_reg);
        % clear data_gpu_reg outs_gpu data_gpu
        % gpuDevice([]);
        
        avg_img = mean(data_g_reg,3);
        f1=figure;
        imagesc(avg_img);
        f1.WindowState = "maximized";
        title(['Select ROI ' expt(day).mouse ' ' pre_or_post]);
        mask_poly = roipoly();
        close();

        roiPixels = avg_img;
        roiPixels(~mask_poly) = 0;
        f2 = figure;
        f2.Name = ['PolyROI_' expt(day).mouse '_' pre_or_post];
        imagesc(roiPixels);
        title(['Polygon ROI ' expt(day).mouse ' ' pre_or_post]);
        
        TC = stackGetTimeCourses(data_g_reg,mask_poly);
        TCmat{iday,1} = TC;
    end
    polyTCcell{id,1} = TCmat;
end

TCs_SST = polyTCcell;

save(fullfile(fn_out_root,'TCs_SST'),'TCs_SST','-v7.3');


%% 
% load(fullfile(fn_out_root,'TCs_SST'));
for iMouse = 1:size(TCs_SST,1)
    this_mouse_data = TCs_SST{iMouse,1};
    mouse_id = TCs_SST{iMouse,2};
    for row = 1:length(this_mouse_data)
        this_mouse_data{row,1} = this_mouse_data{row,1} - mean(this_mouse_data{row,1});
    end

    boCC_day_data = cell2mat(this_mouse_data);
    ymax = max(boCC_day_data);
    ymin = min(boCC_day_data);

    for iday = 1:length(this_mouse_data)
        CCis_TC = this_mouse_data{iday};
        if iday == 1
            pre_or_post = 'pre';
        elseif iday == 2
            pre_or_post = 'post';
        end

        f1 = figure;
        plot(CCis_TC)
        f1.Name = ['PolyTC_' mouse_id '_' pre_or_post];
        sgtitle(['Polygon TC ' mouse_id ' ' pre_or_post]);
        ylim([ymin ymax])
        ylabel('F');
        xlabel('Frame');
    end
end


save_all_open_figs(fullfile(fn_out_root,'plots'));