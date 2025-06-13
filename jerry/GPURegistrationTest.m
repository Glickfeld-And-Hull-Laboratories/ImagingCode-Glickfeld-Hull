clear all; clear global; close all;
clc

dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories


fn_analysis_root = 'G:\home\ACh\Analysis\2p_analysis';
fn_data_root = 'G:\home\ACh\Data\2p_data';
fn_out_root = 'G:\home\jerry\analysis\twophoton\regTest';

session_id_TH = [28]; % enter post-DART session IDs

%% load data and reg info TH


ds = 'DART_expt_info';
eval(ds);
iday = 1;
id = 1;
% for id = 1:length(session_id_TH)
    pre_day = expt(session_id_TH(id)).multiday_matchdays;
    post_day = session_id_TH(id);
    
    days = [pre_day post_day];
    
    % for iday = 1:length(days)
        % load outs
        td = days(iday);
        if iday == 1
            pre_or_post = 'pre';
        elseif iday == 2
            pre_or_post = 'post';
        end

        fn_ana_full = fullfile(fn_analysis_root,expt(td).exptType,expt(td).mouse,expt(td).date,expt(td).contrastxori_runs);
        load(fullfile(cell2mat(fn_ana_full),'regOuts&Img.mat'));
        % find original data
        fn_data_full = fullfile(fn_data_root,expt(td).mouse,expt(td).date,expt(td).contrastxori_runs);
        cd(cell2mat(fn_data_full));
        sesh2load = [cell2mat(expt(td).contrastxori_runs) '_000_000'];
        data_g = sbxread(sesh2load,0,86400);
        data_g = squeeze(data_g(1,:,:,:));
        clear global;

        %first strat
        tic
        [~,data_g_reg] = stackRegister_MA(data_g,[],[],outs);
        t1 = toc;
        % outs_gpu = gpuArray(outs(1:2400,:));
        % data_gpu = gpuArray(data_g);
        % [~,data_gpu_reg] = stackRegister_TH(data_gpu,[],[],outs_gpu);

        % disp(cell2mat(['data is from ' expt(td).mouse ' ' expt(td).date ' ' expt(td).contrastxori_runs]));
        % data_g_reg = gather(data_gpu_reg);
        % clear data_gpu_reg outs_gpu data_gpu
        % gpuDevice([]);
        
        %second strat
        tic
        [~,data_g_reg_gpu] = stackRegister_TH(data_g,[],[],outs);
        t2 = toc;

        

        % full reg
        tic
        [outs_fresh1,data_g_reg] = stackRegister(data_g,regImg);
        t3 = toc;
        figure();
        imagesc(mean(data_g_reg,3))
        sgtitle(['CPU reg fov, t = ' num2str(round(t3,2))]);

        tic
        [outs_fresh2,data_g_reg_gpu] = stackRegister_TH_byFrame(data_g,regImg);
        t4 = toc;
        figure();
        imagesc(mean(data_g_reg_gpu,3))
        sgtitle(['GPU reg fov, t = ' num2str(round(t4,2))])
        
    % end

% end

%% don't run this chunk 
% frames = zeros(550,770,5000);
% allframes = zeros(550,770,90000);
% gpuFrames = zeros(550,770,90000);
% 
% tic
% for i = 1:900
%     t = gpuArray(allframes(:,:,i));
% end
% t1 = toc;
% % Faster method:
% tic
% for j = 1:3
%     temp = gpuArray(frames);
% end
% t2 = toc;
% 
% gpuDevice([]);
% clear all;