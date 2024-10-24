clear all
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);

% redo_exp.date = {'110920' '110817' '110820' '110819' '110817'};
% redo_exp.mouse = {'DR9' 'AC44' 'AC45' 'Y26' 'AC44'};
% redo_exp.userun = {[1:4] [4:6] [1:4] [3:5] [1:3]};
% redo_exp.zoom = {1.5 1 1 1.5 1};
% redo_exp.count_prot = {1 2 1 2 1};
% redo_exp.dirs = {2 2 2 2 2};
% redo_exp.run = {0 0 0 0 0};
% redo_exp.blanks = {1 1 1 1 1};
% redo_exp.area = {'PM' 'AL' 'LM' 'LM' 'LM'};

redo_exp.date = {'111125' '110916'};
redo_exp.mouse = {'Y28' 'Y27'};
redo_exp.userun = {[1:2] [3:5]};
redo_exp.zoom = {1 1.5};
redo_exp.count_prot = {1 2};
redo_exp.dirs = {1 2};
redo_exp.run = {0 0};
redo_exp.blanks = {1 1};
redo_exp.area = {'PM' 'V1'};

for iexp = 1:5
    mouse = char(redo_exp.mouse{iexp});
    date = char(redo_exp.date{iexp});
    userun = redo_exp.userun{iexp};
    count_protocol = redo_exp.count_prot{iexp};
    run = redo_exp.dirs{iexp};
    blanks = redo_exp.blanks{iexp};
    dirs = redo_exp.dirs{iexp};
    zoom = redo_exp.zoom{iexp};
    area = redo_exp.area{iexp};
    
    if area == 'PM'
        iArea = 1;
    elseif area == 'LM'
        iArea = 2;
    elseif area == 'AL'
        iArea = 3;
    end
    
    if dirs ==1
        nCond = 25;
    elseif dirs ==2
        nCond = 50;
    end
   
        
    base = 'G:\users\lindsey\analysisLG\active mice';    
    outDir = fullfile(base, mouse,date);

    fn_lbub = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
    load(fn_lbub);
    fn_fit  = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct.mat']);
    load(fn_fit); 
    i = [];
    j = [];
    fn_local  = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_local_max.mat']);
    load(fn_local); 
    fn_resp = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_resp.mat']);
    load(fn_resp);
    fn_reps= fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
    load(fn_reps);

    nReps = sum(stim_reps(1,:));
    if dirs == 2
        nCond = 50;
        add = 1;
        stim_reps_dir = zeros(1,26);
        for iCond = 1:nCond/2
        stim_reps_dir(1,iCond) = sum(stim_reps(1,add:add+1));
        add = add+2;
        end
        stim_reps_dir(1,end) = stim_reps(1,end);
        stim_reps = stim_reps_dir;
    end

    dF_mat = zeros(size(resp_dFoverF,1), 25);
    start = 1;
    for iCond = 1:25
        nRep = stim_reps(iCond);
        dF_mat(:,iCond) = mean(resp_dF(:,start:start-1+nRep),2);
        start = start+nRep;
    end

    all_fits(iArea).expt(iexp).n = [size(lbub_fits,1) length(goodfit_ind)];


    for iCell = 1:size(lbub_fits,1)
        all_fits(iArea).expt(iexp).bouton(iCell).dF_fit = Fit_struct(iCell).True.s_.x(:,1);
        all_fits(iArea).expt(iexp).bouton(iCell).sigma_SF = Fit_struct(iCell).True.s_.x(:,2);
        all_fits(iArea).expt(iexp).bouton(iCell).sigma_TF = Fit_struct(iCell).True.s_.x(:,3);
        all_fits(iArea).expt(iexp).bouton(iCell).SF_fit = 2.^Fit_struct(iCell).True.s_.x(:,4);
        all_fits(iArea).expt(iexp).bouton(iCell).TF_fit = 2.^Fit_struct(iCell).True.s_.x(:,5);
        all_fits(iArea).expt(iexp).bouton(iCell).xi_fit = Fit_struct(iCell).True.s_.x(:,6);
        all_fits(iArea).expt(iexp).bouton(iCell).speed = 2.^lbub_fits(iCell,5,4)./2.^lbub_fits(iCell,4,4);
        all_fits(iArea).expt(iexp).bouton(iCell).pos = [i(iCell, :) j(iCell, :)];
        all_fits(iArea).expt(iexp).bouton(iCell).plotfit = Fit_struct(iCell).True.s_.k2b_plot;
        all_fits(iArea).expt(iexp).bouton(iCell).dFoverF = Fit_struct(iCell).True.s_.orig;
        all_fits(iArea).expt(iexp).bouton(iCell).dF = dF_mat';
        if find(goodfit_ind == iCell)>0
            all_fits(iArea).expt(iexp).bouton(iCell).goodfit = 1;
        else
            all_fits(iArea).expt(iexp).bouton(iCell).goodfit = 0;
        end
    end
end

for iArea = 1:3
    nexp = all_fits(iArea).nexp;
    all_boutons = [0 0];
    for iexp = 1:nexp
        all_boutons = all_boutons+all_fits(iArea).expt(iexp).n;
    end
    all_fits(iArea).n = all_boutons;
end
        
    fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
    save(fn_out, 'all_fits');
        
        