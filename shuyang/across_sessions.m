%% assign pathnames and datasets to be analyzed/written for moving dot experiments
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'180414_img1005_1','180414_img1007_1','180414_img1008_1','180417_img1008_1',...
    '180419_img1008_1','180423_img1010_1'}; 
days = {'1005-180414_1','1007-180414_1','1008-180414_1','1008-180417_1','1008-180419_1','1010-180423_1'};
image_dest_base    = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
color_code = {'b','r','k'};

%% SECTION I - draw df/f vs. speed across sessions 
dfOvF_spd_all = {}; isp_all = {};
for i = 1: length(sessions)
    image_dest = [image_dest_base sessions{i} '\' sessions{i}];
    img_anal = load([image_dest '_imgAnalysis.mat' ]);
    isp = img_anal.isp;
    dfOvF_spdmean = img_anal.dfOvF_spdmean;
    dfOvF_spd_all = cat(2,dfOvF_spd_all, dfOvF_spdmean);
    isp_all = cat(2,isp_all, isp);
end

%generate matrix for speeds
isp_temp = cell2mat(cellfun(@size,isp_all, 'UniformOutput',0));
isp_length = max(isp_temp);
isp_mat = nan(size(isp_all,2), isp_length);
for n = 1: size(isp_mat,1)
    temp = isp_all{n};
    isp_mat(n,1:size(temp,2)) = temp;
end

%generate matrix for df/f_speedmean
dfOvF_spd_temp = cell2mat(cellfun(@size,dfOvF_spd_all, 'UniformOutput',0));
dfOvF_spd_length = max(dfOvF_spd_temp);
dfOvF_spd_mat = nan(size(dfOvF_spd_all,2), dfOvF_spd_length);
for n = 1: size(dfOvF_spd_mat,1)
    temp = dfOvF_spd_all{n};
    dfOvF_spd_mat(n,1:size(temp,2)) = temp;
end


%% SECTION II - df/f for stay&run across sessions 


%% SECTION III - df/f and speed 300ms across sessions


%% SECTION IV - run trig ave across sessions


%% load things for reverse moving dot experiments
sessions = {'180417_img1005_1','180419_img1005_1','180419_img1007_1','180423_img1005_1',...
    '180424_img1008_1','180425_img1008_1','180428_img1008_1','180429_img1008_1','180430_img1005_1',...
    '180430_img1007_1','180430_img1008_1','180430_img1010_1','180505_img1007_1','180505_img1008_1','180505_img1010_1'};

days = {'1005-180417_1','1005-180419_1','1007-180419_1','1005-180423_1','1008-180424_1','1008-180425_1','1008-180428_1',...
    '1008-180429_1','1005-180430_1','1007-180430_1','1008-180430_1','1010-180430_1','1007-180505_1','1008-180505_1',...
    '1010-180505_1'};


%% SECTION V - reverse trig ave across sessions


%% SECTION VI - reverse trig ave during stay across sessions


%% SECTION VII reverse trig ave during run across sessions