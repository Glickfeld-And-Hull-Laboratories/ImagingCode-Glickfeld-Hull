
clear all; clear global; close all
clc
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};

eval(ds);
% 259 269 emx atropine - forward matched
% 206 210 214 atropine SOM
% 201 197 emx YM90K - forward matched
% 178 190 294 %good quality SOM YM90K
%138 142 163 171 178 190 294 307 for retreat talk
%294 307 323 NES with DART
%peg 303 311 319 329
%Oct 2023: 138 142 163 171 178 190 294 307 323 333


sess_list = [138 142 163 171 178 190 294 307 333 323 303 311 319 329 355 359];%enter all the sessions you want to concatenate
nSess=length(sess_list);

nd=2;%hard coding for two days per experimental session

% INDICATE THE PRE VS. POST DAYS DEPENDING ON THE ORDER OF MATCHING
prompt = 'Which sesson was used as reference for matching: 0- baseline, 1- post-DART';
            x = input(prompt);
            switch x
                case 0
                    pre=1; %baeline session, used as reference, is in the 1st position
                    post=2;
                    "baseline used as reference"
                case 1
                  pre=2;
                  post=1; %post-DART session, used as reference, is in the 1st position  
                  "post-DART used as reference"
            end
clear x prompt

targetCon = [.25 .5 1]%what contrast to extract for all data - must be one that all datasets had

frame_rate = 15;

sess_title = string(sess_list(1));
for iSess = 2:nSess
    sess_title = strcat(sess_title,'_',string(sess_list(iSess)));
end
d=string(datetime('today'));

if computer == 'GLNXA64'
    isilonName =  '/home/cc735@dhe.duke.edu/GlickfeldLabShare';
    base = fullfile(isilonName, '/All_Staff/home/ACh/Data/2p_data');
    
else
    isilonName = 'Z:';
    base = fullfile(isilonName, '\home\ACh\Analysis\2p_analysis');
      
end


if nSess == 1
         if expt(sess_list(1)).multiday_timesincedrug_hours>0
            dart_str = [expt(sess_list(1)).drug '_' num2str(expt(sess_list(1)).multiday_timesincedrug_hours) 'Hr'];
        else
            dart_str = 'control';
        end
        
        fnout = fullfile(base,expt(sess_list(1)).mouse,['multiday_' dart_str],d);
else
    fnout= fullfile(base,strcat('concat', sess_title),d);
end
mkdir(fnout);
cd(fnout)
clear d sess_title

zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
%% concatenating data
nCon = length(targetCon)


dataTableConat=[];
nKeep_concat=[];
red_concat = [];
green_concat = [];

cellID_adjustment=0;
for iSess = 1:nSess
    day_id = sess_list(iSess)
    mouse = expt(day_id).mouse;
    thisDrug = expt(day_id).drug;
    drug{iSess}=thisDrug;

    if iSess > 1
        cellID_adjustment=max(temp_table.cell_ID_unique); %this should get saved until the next loop;
    end

    if expt(day_id).multiday_timesincedrug_hours>0
        dart_str = [expt(day_id).drug '_' num2str(expt(day_id).multiday_timesincedrug_hours) 'Hr'];
    else
        dart_str = 'control';
    end
    fn_multi = fullfile(base,mouse,['multiday_' dart_str]);


   temp_table =readtable(fullfile(fn_multi,'dataTable.csv'));
    
    temp_table.z_speed=zscor_xnan(temp_table.speed);
    temp_table.z_pupil=zscor_xnan(temp_table.pupil);
    temp_table.cell_ID_unique=temp_table.cellID + cellID_adjustment;

    dataTableConat=[dataTableConat; temp_table];

    load(fullfile(fn_multi,'tc_keep.mat'));
    nKeep = size(tc_trial_avrg_stat{post},2);
    nKeep_concat = [nKeep_concat,nKeep];
    red_concat = [red_concat, red_keep_logical];
    green_concat = [green_concat, green_keep_logical];




end
%
nKeep_total = sum(nKeep_concat);
red_ind_concat = find(red_concat);
green_ind_concat = find(green_concat);

%% mixed model - making dataset

dataTableConat.day=categorical(dataTableConat.day);
dataTableConat.mouseID=categorical(dataTableConat.mouseID);
dataTableConat.cellID=categorical(dataTableConat.cell_ID_unique);
dataTableConat.cellType=categorical(dataTableConat.cellType);
dataTableConat.drug = categorical(dataTableConat.drug);
dataTableConat.z_speed=zscor_xnan(dataTableConat.speed);
dataTableConat.z_pupil=zscor_xnan(dataTableConat.pupil);


%% mixed model on whole population

%lme = fitlme(dataTableConat,'dfof~cellType+contrast+speed+pupil+day+drug+(cellType*day*drug)+(1|mouseID)+(1|cell_ID_unique:mouseID)')
lme = fitlme(dataTableConat,['dfof~cellType+contrast+z_speed+z_pupil+day+drug+(cellType*day*drug)+' ...
    '(1|mouseID)+(1|cell_ID_unique:mouseID)+(1|day:cell_ID_unique)+(1|drug:cell_ID_unique)'])

%% post-hoc comaparisons for mixed model
%what is the mean df/f for each cell type on each day
vars = ["dfof"];
factors = ["drug","cellType","day"];
meanScoresByFactor = varfun(@mean, ...
                            dataTableConat, ...
                            "InputVariables",vars, ...
                            "GroupingVariables",factors)

%% mixed model cell by cell
betas = nan(nKeep_total,5);
for iCell = 1:nKeep_total
    inds = find(dataTableConat.cell_ID_unique==iCell);
    temp_data = dataTableConat(inds,:);
    lm_cell = fitlm(temp_data,'dfof~direction+contrast+z_speed+z_pupil+day');
    betas(iCell,:)=lm_cell.Coefficients.Estimate(2:6);
    %betas(iCell,5)=betas(iCell,5)*-1;
end

%
%betas: 1=direction, 2=contrast, 3=z-score speed, 4=z-score pupil, 5=day (pre/post),
%%
%make intersectional indices of cell type x drug
DART_cells = unique(dataTableConat.cell_ID_unique(find(dataTableConat.drug=='YM90K-DART'))); %find the unique cell IDs in the rows of the table where drug is DART
PEG_cells =  unique(dataTableConat.cell_ID_unique(find(dataTableConat.drug=='YM90K-PEG'))); %likewise for PEG

DART_green = intersect(green_ind_concat,DART_cells);
DART_red =  intersect(red_ind_concat,DART_cells);

PEG_green = intersect(green_ind_concat,PEG_cells);
PEG_red =  intersect(red_ind_concat,PEG_cells);


figure; 
subplot(2,1,1)
histogram(betas(DART_green,5));hold on;histogram(betas(DART_red,5));title('betas fore pre/post'); title('DART');xlim([-.15 .15]);xlabel('beta for day');hold off
subplot(2,1,2)
histogram(betas(PEG_green,5));hold on;histogram(betas(PEG_red,5));title('betas fore pre/post'); title('PEG');xlim([-.15 .15]);xlabel('beta for day');hold off