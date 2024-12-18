clear all

%% expt params
date = '111218';
mouse = 'M32';
userun = [1:5];

count_protocol = 1;
blanks = 1;
run = 0;

nCond =25;
%%
P = 1;
% nON = 10;
% nOFF = 10;
% nPlanes = 1;
% begin = 7;
% TFSFetc = [1:2];
% pre_win = [1 4];
% post_win = [5 14];
nON = 10;
nOFF = 40;
nPlanes = 1;
begin = 31;
TFSFetc = [1:2];
pre_win = [1 10];
post_win = [11 20];

base = 'G:\users\lindsey\analysisLG\active mice';
running_base = 'G:\users\lindsey\dataLG\Running data';
outDir = fullfile(base, mouse,date);
%% load image for choosing rois
fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_max_blur.tif']);
max_dF = readtiff(fn);
bwimgarea = imCellEditInteractive2(max_dF,[]);
area_mask = bwlabel(bwimgarea);
figure; imagesc(area_mask);

a = 'LM'; b = 'AL'; c = 'RL'; d = 'V1';  e = 'A'; f = 'PM'; g = 'AM';
area_list = strvcat(a,b,c,d,e,f,g);

fn_out = fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_area_mask.mat']);
save(fn_out, 'area_mask', 'area_list');

%% resort expt
%load params and creat seq file
eval(['PARAMS_' date '_' mouse]);
resort_seq_only

%sort and average stacks
if blanks ==0;
    stack_sort_noblanks
else
    stack_sort
end;
%measure and sort by running
Running_LG

%% clean up images
%normalize luminance
Epi_clean

%% get timecourses

roi_sort_norm

%% plot data
%timecourses
area_colors
       
Make_TC_subplots

%% stats
fn_area = fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_area_mask.mat']);
fn_var = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_vars.mat']);
fn_reps = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
fn_reps_run = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stim_reps_run.mat']);
fn_reps_norun = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stim_reps_norun.mat']);
fn_resp = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp_POST' num2str(post_win) '.mat']);
fn_resp_run = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp_POST' num2str(post_win) '_run_norun.mat']);
fn_resp_norm = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp_norm_POST' num2str(post_win) '.mat']);
fn_avg = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_roi_avg.mat']);
load(fn_area);
load(fn_var);
load(fn_reps);
load(fn_resp);
load(fn_resp_norm);
load(fn_avg)
if run == 1;
    load(fn_reps_run);
    load(fn_reps_norun);
    load(fn_resp_run);
end

areas = size(area_list,1);
Nshuf = 0; %0 means only run for true data
Cond_USE = 1:nCond;
nCells = areas;
SFTF_fit_LG
