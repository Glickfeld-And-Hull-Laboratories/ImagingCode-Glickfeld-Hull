clear all;
close all;
clc;
doRedChannel = 0;
ds = 'F0F1F2_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
fnout = fullfile(LG_base,'Analysis','2P');
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'PhaseRev');
svName = 'F0F1F2';
dirF0_all = [];
dirF1_all = [];
prF1_all = [];
prF2_all = [];
mouse_list = [];
    
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    mouse_list = strvcat(mouse_list, mouse);
    date = expt(iexp).date;
    datemouse = [date '_' mouse];
    dirFolder = expt(iexp).dirFolder;
    nrun = length(dirFolder);
    dir_run_str = catRunName(cell2mat(dirFolder), nrun);
    datemousedirrun = [date '_' mouse '_' dir_run_str];
    prFolder = expt(iexp).prFolder;
    nrun = length(prFolder);
    pr_run_str = catRunName(cell2mat(prFolder), nrun);
    datemouseprrun = [date '_' mouse '_' pr_run_str];
    fprintf([mouse ' ' date '\n'])
    
    dirdata = load(fullfile(fnout, datemouse, datemousedirrun, [datemousedirrun '_dirAnalysis.mat']));
    prdata = load(fullfile(fnout, datemouse, datemouseprrun, [datemouseprrun '_f1f2.mat']));
    
    dirF0_all = [dirF0_all dirdata.f0];
    dirF1_all = [dirF1_all dirdata.f1];
    prF1_all = [prF1_all prdata.f1];
    prF2_all = [prF2_all prdata.f2];
end
 
resp_ind_all = find(dirF0_all>0.02 & prF1_all>0.02);

figure;
subplot(2,2,1)
scatter(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all),prF2_all(resp_ind_all)./prF1_all(resp_ind_all))
xlabel('F1/F0')
ylabel('F2/F1')
xlim([0 0.5])
ylim([0 2])
subplot(2,2,2)
scatter(dirF1_all(resp_ind_all),prF1_all(resp_ind_all))
xlabel('drifting F1')
ylabel('contrast rev F1')
xlim([0 1])
ylim([0 1])
refline(1)
subplot(2,2,3)
histogram(dirF1_all(resp_ind_all)./dirF0_all(resp_ind_all),[0:0.1:2])
xlabel('F1/F0')
xlim([0 2])
subplot(2,2,4)
histogram(prF2_all(resp_ind_all)./prF1_all(resp_ind_all),[0:0.1:2])
xlabel('F2/F1')
xlim([0 2])
sgtitle(['n = ' num2str(length(resp_ind_all)) ' cells; ' num2str(size(mouse_list,1)) ' mice'])
print(fullfile(summaryDir,[svName '_summary.pdf']),'-dpdf','-fillpage')
