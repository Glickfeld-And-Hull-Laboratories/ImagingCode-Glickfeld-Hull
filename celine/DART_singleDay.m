clear all; clear global;  close all
clc

ds = 'DART_V1_atropine_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);
doGreenOnly = true;
doCorrImg = true;

day_id =64;
experimentFolder = 'SST_atropine';

if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    base = fullfile('/All_Staff/home/ACh/Analysis/2p_analysis',experimentFolder);
    beh_prefix = strcat(isilonName,'/All_Staff/Behavior/Data/');
elseif string(hostname) == 'NB-NUKE'
    isilonName = 'Z:/All_Staff';
    base = fullfile('/home/ACh/Analysis/2p_analysis/',experimentFolder);
    beh_prefix = strcat('Z:/All_Staff/Behavior/Data/');
else
    isilonName = '';
    base = fullfile('/home/ACh/Analysis/2p_analysis',experimentFolder);
    beh_prefix = strcat('Z:\Behavior\Data\');
end
%% load TCs
mouse = expt(day_id).mouse;
expDate = expt(day_id).date;
run = char(eval(['expt(day_id).' cell2mat(dataStructLabels) '_runs']));
fn = fullfile(isilonName,base,mouse,expDate,run);
dat = 'data-';
cd(fn);
load('TCs.mat')
