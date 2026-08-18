%% read output files from ImageJ
clear all; close all
clc

isilonName = 'G:';
base = fullfile(isilonName, '\home\jerry\analysis\ilastik\DARTquant\csvs\');
% metadata-directed file access to be added here
fn = fullfile(base, 'TH015_03_r4_m4.csv');

dat = readtable(fn);

fprintf(['A total of ' num2str(size(dat,1)) ' ROIs were identified before cleanup\n'])

%% clean up ROIs by several criteria
dat_a = table2array(dat);

del = dat_a(:,2) < 10; %find ROIs smaller than 10 pix and delete
dat_a(del,:) = [];

histogram(dat_a(:,2))