clear all; clear global; close all;
clc

dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories


fn_analysis_root = 'G:\home\ACh\Analysis\2p_analysis';
fn_data_root = 'G:\home\ACh\Data\2p_data';
fn_out_root = 'G:\home\ACh\Analysis\2p_analysis\epileptiform_analysis\movies';

session_id_TH = [28 31 34]; % enter post-DART session IDs for movies to be made;
session_id_CC = [169 177 183]; % enter post-DART session IDs for movies to be made;
pre_mov_TH = 1;
post_mov_TH = 1;
pre_mov_CC = 1;
post_mov_CC = 1;

%% movie - TH
ds = 'DART_expt_info';
eval(ds);

startFrame = 151;
totFrames = 450;
endFrame = startFrame+totFrames-1;
trimmedStart = startFrame + 90;



for id = 1:length(session_id_TH)

    if pre_mov_TH == 1
        pre_day = expt(session_id_TH(id)).multiday_matchdays;
    else
        pre_day = [];
    end
    
    if post_mov_TH == 1
        post_day = session_id_TH(id);
    else
        post_day = [];
    end
    
    days = [pre_day post_day];
    
    for iday = 1:length(days)
        % load outs
        day = days(iday);
        if length(days) == 2
            if iday == 1
                pre_or_post = 'pre';
            elseif iday == 2
                pre_or_post = 'post';
            end
        elseif length(days) == 1
            if post_mov_TH == 1
                pre_or_post = 'post';
            elseif pre_mov_TH == 1
                pre_or_post = 'pre';
            end
        end

        fn_ana_full = fullfile(fn_analysis_root,expt(day).exptType,expt(day).mouse,expt(day).date,expt(day).contrastxori_runs);
        load(fullfile(cell2mat(fn_ana_full),'regOuts&Img.mat'));
        % find original data
        fn_data_full = fullfile(fn_data_root,expt(day).mouse,expt(day).date,expt(day).contrastxori_runs);
        cd(cell2mat(fn_data_full));
        sesh2load = [cell2mat(expt(day).contrastxori_runs) '_000_000'];
        data_g = sbxread(sesh2load,startFrame-1,totFrames);
        data_g = squeeze(data_g(1,:,:,:));
        clear global
        
        outs_gpu = gpuArray(outs(startFrame:endFrame,:));
        data_gpu = gpuArray(data_g);
        [~,data_gpu_reg] = stackRegister_TH(data_gpu,[],[],outs_gpu);
        disp(['data is from' expt(day).mouse ' ' expt(day).date]);
        max_data = max(data_gpu_reg,[],3);
        max_data_cpu = gather(max_data);
        sorted = sort(max_data_cpu(:));
        [N,edges] = histcounts(max_data_cpu,20);
        data_g_reg = gather(data_gpu_reg);
        clear data_gpu_reg outs_gpu max_data data_gpu
        gpuDevice([]);
        
        figure;
        avg_img = mean(data_g_reg,3);
        imagesc(avg_img)
        
        normA = double(data_g_reg);
        normA = normA - min(normA(:));
        normA = normA ./ double(sorted(end-1000));
        normA(normA>1) = 1;
        
        % ---- marker parameters ----
        N_marker = 90;     % Show marker every 90 frames
        duration = 30;      % Marker persists for 30 frames
        markerPos = [60, 60];            % x, y
        radius = 20;                      % circle radius
        
        % ---- Write video with marker ----
        movie_out = fullfile(fn_out_root,[expt(day).mouse pre_or_post '_movie']);
        v = VideoWriter(movie_out, 'MPEG-4');
        v.FrameRate = 15;
        open(v);
        
        for k = 1:size(normA, 3)
            % Convert frame to uint8 RGB
            frame_gray = uint8(normA(31:end,31:end,k) * 255);
            frame_rgb = repmat(frame_gray, [1 1 3]);
            % frame_rgb_down = imresize(frame_rgb,[500 766]);
            % frame_rgb_gauss = imgaussfilt(frame_rgb,0.5);
        
            % Add marker every N frames
            if mod(k, N_marker) < duration
                frame_rgb = insertShape(frame_rgb, 'FilledCircle', [markerPos, radius], ...
                                        'Color', 'red', 'Opacity', 1);
            end
            writeVideo(v, frame_rgb);
        end
        close(v);

        inputFile = [movie_out '.mp4'];
        outputFile = [movie_out '_trimmed_' 's' num2str(trimmedStart) '_e' num2str(endFrame)];

        vidReader = VideoReader(inputFile);
        frameRate = vidReader.FrameRate;

        vidWriter = VideoWriter(outputFile, 'MPEG-4');
        vidWriter.FrameRate = frameRate;
        open(vidWriter);

        startTime = 6; % start after glitch
        endTime = 36;

        while hasFrame(vidReader)
        frameTime = vidReader.CurrentTime;
    
        if frameTime >= startTime && frameTime <= endTime
            frame = readFrame(vidReader);
            writeVideo(vidWriter, frame);
        elseif frameTime > endTime
            break;
        else
            % Still need to read frame to move forward in time
            readFrame(vidReader); 
        end
        end

        close(vidWriter);
        % close(v);
    end

end
%% movie - CC
ds = 'DART_V1_contrast_ori_Celine';
eval(ds);

startFrame = 151;
totFrames = 450;
endFrame = startFrame+totFrames-1;
trimmedStart = startFrame + 90;



for id = 1:length(session_id_CC)

    if pre_mov_CC == 1
        pre_day = expt(session_id_CC(id)).multiday_matchdays;
    else
        pre_day = [];
    end
    
    if post_mov_CC == 1
        post_day = session_id_CC(id);
    else
        post_day = [];
    end
    
    days = [pre_day post_day];
    
    for iday = 1:length(days)
        % load outs
        day = days(iday);
        if length(days) == 2
            if iday == 1
                pre_or_post = 'pre';
            elseif iday == 2
                pre_or_post = 'post';
            end
        elseif length(days) == 1
            if post_mov_CC == 1
                pre_or_post = 'post';
            elseif pre_mov_CC == 1
                pre_or_post = 'pre';
            end
        end

        fn_ana_full = fullfile(fn_analysis_root,expt(day).exptType,expt(day).mouse,expt(day).date,expt(day).contrastxori_runs);
        load(fullfile(cell2mat(fn_ana_full),'regOuts&Img.mat'));
        % find original data
        fn_data_full = fullfile(fn_data_root,expt(day).mouse,expt(day).date,expt(day).contrastxori_runs);
        cd(cell2mat(fn_data_full));
        sesh2load = [cell2mat(expt(day).contrastxori_runs) '_000_000'];
        data_g = sbxread(sesh2load,startFrame-1,totFrames);
        data_g = squeeze(data_g(1,:,:,:));
        clear global
        
        outs_gpu = gpuArray(outs(startFrame:endFrame,:));
        data_gpu = gpuArray(data_g);
        [~,data_gpu_reg] = stackRegister_TH(data_gpu,[],[],outs_gpu);
        disp(['data is from' expt(day).mouse ' ' expt(day).date]);
        max_data = max(data_gpu_reg,[],3);
        max_data_cpu = gather(max_data);
        sorted = sort(max_data_cpu(:));
        [N,edges] = histcounts(max_data_cpu,20);
        data_g_reg = gather(data_gpu_reg);
        clear data_gpu_reg outs_gpu max_data data_gpu
        gpuDevice([]);
        
        figure;
        avg_img = mean(data_g_reg,3);
        imagesc(avg_img)
        
        normA = double(data_g_reg);
        normA = normA - min(normA(:));
        normA = normA ./ double(sorted(end-1000));
        normA(normA>1) = 1;
        
        % ---- marker parameters ----
        N_marker = 90;     % Show marker every 90 frames
        duration = 30;      % Marker persists for 30 frames
        markerPos = [60, 60];            % x, y
        radius = 20;                      % circle radius
        
        % ---- Write video with marker ----
        movie_out = fullfile(fn_out_root,[expt(day).mouse pre_or_post '_movie']);
        v = VideoWriter(movie_out, 'MPEG-4');
        v.FrameRate = 15;
        open(v);
        
        for k = 1:size(normA, 3)
            % Convert frame to uint8 RGB
            frame_gray = uint8(normA(31:end,31:end,k) * 255);
            frame_rgb = repmat(frame_gray, [1 1 3]);
            % frame_rgb_down = imresize(frame_rgb,[500 766]);
            % frame_rgb_gauss = imgaussfilt(frame_rgb,0.5);
        
            % Add marker every N frames
            if mod(k, N_marker) < duration
                frame_rgb = insertShape(frame_rgb, 'FilledCircle', [markerPos, radius], ...
                                        'Color', 'red', 'Opacity', 1);
            end
            writeVideo(v, frame_rgb);
        end
        close(v);

        inputFile = [movie_out '.mp4'];
        outputFile = [movie_out '_trimmed_' 's' num2str(trimmedStart) '_e' num2str(endFrame)];

        vidReader = VideoReader(inputFile);
        frameRate = vidReader.FrameRate;

        vidWriter = VideoWriter(outputFile, 'MPEG-4');
        vidWriter.FrameRate = frameRate;
        open(vidWriter);

        startTime = 6; % start after glitch
        endTime = 36;

        while hasFrame(vidReader)
        frameTime = vidReader.CurrentTime;
    
        if frameTime >= startTime && frameTime <= endTime
            frame = readFrame(vidReader);
            writeVideo(vidWriter, frame);
        elseif frameTime > endTime
            break;
        else
            % Still need to read frame to move forward in time
            readFrame(vidReader); 
        end
        end

        close(vidWriter);
        % close(v);
    end

end