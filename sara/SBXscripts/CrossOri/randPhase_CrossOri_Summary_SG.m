close all; clear all; clc;
doRedChannel = 0;
SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');

doPlot = 1;
ds = ['CrossOriRandPhase_15Hz_ExptList_SG'];
svName = 'randPhase';
eval(ds)
driver = 'SLC';
img_area = {'LM';'L2/3'}; %LM
inj_area = 'LM';

rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);

sigCells = [];
fractSigCells = [];
respCells = [];
resp_ind_all = [];

p_anova_all_all = [];
b_all_all = [];
amp_all_all = [];
pha_all_all = [];
per_all_all = [];
sse_all_all = [];
R_square_all_all = [];
yfit_all_all = [];

p_anova_shuf_all = [];
b_shuf_all = [];
amp_shuf_all = [];
pha_shuf_all = [];
per_shuf_all = [];
sse_shuf_all = [];
R_square_shuf_all = [];
yfit_shuf_all = [];

p_anova_downsamp_all = [];
b_downsamp_all = [];
amp_downsamp_all = [];
pha_downsamp_all = [];
per_downsamp_all = [];
sse_downsamp_all = [];
R_square_downsamp_all = [];
yfit_downsamp_all = [];

p_anova_thresh_all = [];
b_thresh_all = [];
amp_thresh_all = [];
pha_thresh_all = [];
per_thresh_all = [];
R_square_thresh_all = [];
sse_thresh_all = [];

p_anova_shuf_thresh_all = [];
b_shuf_thresh_all = [];
amp_shuf_thresh_all = [];
pha_shuf_thresh_all = [];
per_shuf_thresh_all = [];
R_square_shuf_thresh_all = [];
sse_shuf_thresh_all = [];

test_thresh_resp_all = [];
mask_thresh_resp_all = [];

maskPhas_all_all = [];

plaidSI_all = [];
testPI_all = [];
plaidSI_thresh_all = [];
testPI_thresh_all = [];
OSI_all = [];
max_dir_all = [];

resp_avg_max_all = [];
test_avg_all = [];
mask_avg_all = [];


mouse_list = [];

totCells = zeros(nexp,1);

for iexp = [47]
%     if strcmp(expt(iexp).inj_loc,inj_area) & strcmp(expt(iexp).driver,driver)        
        mouse = expt(iexp).mouse;
        mouse_list = strvcat(mouse_list, mouse);
        date = expt(iexp).date;
        
        if isfield(expt,'copFolder') 
            ImgFolder = expt(iexp).copFolder;
        else
            ImgFolder = expt(iexp).coFolder;
        end
        nrun = length(ImgFolder);
%         run_str = catRunName(cell2mat(ImgFolder), nrun);
        run_str = 'runs-002';

        if strcmp(expt(iexp).saveLoc,'sara')
            SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
        elseif strcmp(expt(iexp).saveLoc,'lindsey')
            SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
        end
        
        fprintf([mouse ' ' date '\n'])
        
        
        load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
        load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
    

        nCells = size(resp_cell{end,end,end},1);
        totCells(iexp,:) = nCells;

        p_anova_all(find(p_anova_all==0)) = NaN;
        p_anova_all_all = [p_anova_all_all;  p_anova_all];
        b_all_all = [b_all_all; b_hat_all];
        amp_all_all = [amp_all_all; amp_hat_all];
        pha_all_all = [pha_all_all; pha_hat_all];
        per_all_all = [per_all_all; per_hat_all];
        R_square_all_all = [R_square_all_all; R_square_all];
        sse_all_all = [sse_all_all; sse_all];

        p_anova_shuf(find(p_anova_shuf==0)) = NaN;
        p_anova_shuf_all = [p_anova_shuf_all;  p_anova_shuf];
        b_shuf_all = [b_shuf_all; b_hat_shuf];
        amp_shuf_all = [amp_shuf_all; amp_hat_shuf];
        pha_shuf_all = [pha_shuf_all; pha_hat_shuf];
        per_shuf_all = [per_shuf_all; per_hat_shuf];
        R_square_shuf_all = [R_square_shuf_all; R_square_shuf];
        sse_shuf_all = [sse_shuf_all; sse_shuf];
       
        
        yfit_all_all = cat(1,yfit_all_all, yfit_all);
        yfit_shuf_all = cat(1,yfit_shuf_all, yfit_shuf);
        
        yfit_all_all_max = max(yfit_all_all,[],2);
        amp_all_all_norm = (amp_all_all ./ yfit_all_all_max);
        b_all_all_norm = (b_all_all ./ yfit_all_all_max);


        resp_ind_all = [resp_ind_all; resp_ind+sum(totCells(1:iexp-1,:),1)];

        plaid_resp = mean(resp_cell{end,end,1},2);
        mask_resp = mean(resp_cell{end,1,1},2);
        test_resp = mean(resp_cell{1,end,1},2);
        plaid_resp(find(plaid_resp<0)) = 0;
        mask_resp(find(mask_resp<0)) = 0;
        test_resp(find(test_resp<0)) = 0;

        plaidSI_all = [plaidSI_all; (plaid_resp-(mask_resp+test_resp)) ./ (plaid_resp + mask_resp + test_resp)];
        testPI_all = [testPI_all; abs((test_resp-mask_resp) ./ (mask_resp+test_resp))];
        
        plaid_thresh_resp = mean(resp_cell{end,end,1},2).*0.75;
        mask_thresh_resp = mean(resp_cell{end,1,1},2).*0.75;
        test_thresh_resp = mean(resp_cell{1,end,1},2).*0.75;
        plaid_thresh_resp(find(plaid_thresh_resp<0.04)) = 0;
        mask_thresh_resp(find(mask_thresh_resp<0.04)) = 0;
        test_thresh_resp(find(test_thresh_resp<0.04)) = 0;
        
        test_thresh_resp_all = [test_thresh_resp_all; test_thresh_resp];
        mask_thresh_resp_all = [mask_thresh_resp_all; mask_thresh_resp];
        
        plaidSI_thresh_all = [plaidSI_thresh_all; (plaid_thresh_resp-(mask_thresh_resp+test_thresh_resp)) ./ (plaid_thresh_resp + mask_thresh_resp + test_thresh_resp)];
        testPI_thresh_all = [testPI_thresh_all; abs((test_thresh_resp-mask_thresh_resp) ./ (mask_thresh_resp+test_thresh_resp))];
        
        maskPhas_all_all = [maskPhas_all_all; maskPhas_all];
        
        if exist('resp_avg_max', 'var')
            resp_avg_max_all = [resp_avg_max_all; resp_avg_max];
        else
        end
%     end
end
save(fullfile(summaryDir,[svName '_Summary_' inj_area '_' driver '_0' num2str(iexp) '.mat']), 'p_anova_all_all','b_all_all','amp_all_all','R_square_all_all', 'R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','pha_all_all', 'OSI_all', 'max_dir_all','mouse_list', 'maskPhas_all_all')

%% summary for resp plaid

close all; clear all; clc;
doRedChannel = 0;
SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');

doPlot = 1;
ds = ['CrossOriRandPhase_15Hz_ExptList_SG'];
svName = 'randPhase';
eval(ds)
driver = 'SLC';
area = 'V1'; %LM


rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);

sigCells = [];
fractSigCells = [];
respCells = [];
resp_ind_all = [];

p_anova_all_all = [];
b_all_all = [];
amp_all_all = [];
pha_all_all = [];
per_all_all = [];
sse_all_all = [];
R_square_all_all = [];
yfit_all_all = [];

p_anova_shuf_all = [];
b_shuf_all = [];
amp_shuf_all = [];
pha_shuf_all = [];
per_shuf_all = [];
sse_shuf_all = [];
R_square_shuf_all = [];
yfit_shuf_all = [];

p_anova_downsamp_all = [];
b_downsamp_all = [];
amp_downsamp_all = [];
pha_downsamp_all = [];
per_downsamp_all = [];
sse_downsamp_all = [];
R_square_downsamp_all = [];
yfit_downsamp_all = [];

p_anova_thresh_all = [];
b_thresh_all = [];
amp_thresh_all = [];
pha_thresh_all = [];
per_thresh_all = [];
R_square_thresh_all = [];
sse_thresh_all = [];

p_anova_shuf_thresh_all = [];
b_shuf_thresh_all = [];
amp_shuf_thresh_all = [];
pha_shuf_thresh_all = [];
per_shuf_thresh_all = [];
R_square_shuf_thresh_all = [];
sse_shuf_thresh_all = [];

test_thresh_resp_all = [];
mask_thresh_resp_all = [];

maskPhas_all_all = [];

plaid_resp_all = [];
plaid_std_all = [];
plaid_cov_all = [];

plaidSI_all = [];
testPI_all = [];
plaidSI_thresh_all = [];
testPI_thresh_all = [];
OSI_all = [];
max_dir_all = [];

resp_avg_max_all = [];
test_avg_all = [];
mask_avg_all = [];
respplaid_ind_all = [];

respPlaid_avg_all = [];
respTest_max_all = [];
respmask_max_all = [];

mouse_list = [];

totCells = zeros(nexp,1);

for iexp = [8]
        mouse = expt(iexp).mouse;
        mouse_list = strvcat(mouse_list, mouse);
        date = expt(iexp).date;
        
        if isfield(expt,'copFolder') 
            ImgFolder = expt(iexp).copFolder;
        else
            ImgFolder = expt(iexp).coFolder;
        end
        nrun = length(ImgFolder);
        run_str = 'runs-002';;
        
        if strcmp(expt(iexp).saveLoc,'sara')
            SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
        elseif strcmp(expt(iexp).saveLoc,'lindsey')
            SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
        end
        
        fprintf([mouse ' ' date '\n'])
        
        load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
        load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
        load(fullfile(SG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits_respPlaid.mat']))
    

        nCells = size(resp_cell{end,end,end},1);
        totCells(iexp,:) = nCells;

        p_anova_all(find(p_anova_all==0)) = NaN;
        p_anova_all_all = [p_anova_all_all;  p_anova_all];
        b_all_all = [b_all_all; b_hat_all];
        amp_all_all = [amp_all_all; amp_hat_all];
        pha_all_all = [pha_all_all; pha_hat_all];
        per_all_all = [per_all_all; per_hat_all];
        R_square_all_all = [R_square_all_all; R_square_all];
        sse_all_all = [sse_all_all; sse_all];

        p_anova_shuf(find(p_anova_shuf==0)) = NaN;
        p_anova_shuf_all = [p_anova_shuf_all;  p_anova_shuf];
        b_shuf_all = [b_shuf_all; b_hat_shuf];
        amp_shuf_all = [amp_shuf_all; amp_hat_shuf];
        pha_shuf_all = [pha_shuf_all; pha_hat_shuf];
        per_shuf_all = [per_shuf_all; per_hat_shuf];
        R_square_shuf_all = [R_square_shuf_all; R_square_shuf];
        sse_shuf_all = [sse_shuf_all; sse_shuf];
       
        
        yfit_all_all = cat(1,yfit_all_all, yfit_all);
        yfit_shuf_all = cat(1,yfit_shuf_all, yfit_shuf);
        
        yfit_all_all_max = max(yfit_all_all,[],2);
        yfit_all_all_max_rect = yfit_all_all_max;
        yfit_all_all_max_rect(find(yfit_all_all_max<0)) = NaN;
%         amp_all_all_norm = (amp_all_all ./ yfit_all_all_max);
%         b_all_all_norm = (b_all_all ./ yfit_all_all_max);


        resp_ind_all = [resp_ind_all; resp_ind+sum(totCells(1:iexp-1,:),1)];

        plaid_resp = mean(resp_cell{end,end,1},2);
        plaid_resp_std = std(resp_cell{end,end,1},0,2);
        plaid_resp_cov = plaid_resp_std./plaid_resp;
        mask_resp = mean(resp_cell{end,1,1},2);
        test_resp = mean(resp_cell{1,end,1},2);
        plaid_resp(find(plaid_resp<0)) = 0;
        mask_resp(find(mask_resp<0)) = 0;
        test_resp(find(test_resp<0)) = 0;

        plaidSI_all = [plaidSI_all; (plaid_resp-(mask_resp+test_resp)) ./ (plaid_resp + mask_resp + test_resp)];
        testPI_all = [testPI_all; abs((test_resp-mask_resp) ./ (mask_resp+test_resp))];
        
        plaid_thresh_resp = mean(resp_cell{end,end,1},2).*0.75;
        mask_thresh_resp = mean(resp_cell{end,1,1},2).*0.75;
        test_thresh_resp = mean(resp_cell{1,end,1},2).*0.75;
        plaid_thresh_resp(find(plaid_thresh_resp<0.04)) = 0;
        mask_thresh_resp(find(mask_thresh_resp<0.04)) = 0;
        test_thresh_resp(find(test_thresh_resp<0.04)) = 0;
        
        test_thresh_resp_all = [test_thresh_resp_all; test_thresh_resp];
        mask_thresh_resp_all = [mask_thresh_resp_all; mask_thresh_resp];
        
        plaidSI_thresh_all = [plaidSI_thresh_all; (plaid_thresh_resp-(mask_thresh_resp+test_thresh_resp)) ./ (plaid_thresh_resp + mask_thresh_resp + test_thresh_resp)];
        testPI_thresh_all = [testPI_thresh_all; abs((test_thresh_resp-mask_thresh_resp) ./ (mask_thresh_resp+test_thresh_resp))];
        
        maskPhas_all_all = [maskPhas_all_all; maskPhas_all];
        
        resp_avg_max_all = [resp_avg_max_all; resp_avg_max];
        test_avg_all = [test_avg_all; test_avg];
        mask_avg_all = [mask_avg_all; mask_avg];
        respplaid_ind_all = [respplaid_ind_all; respplaid_ind];
        
        plaid_resp_all = [plaid_resp_all; plaid_resp];
        plaid_std_all = [plaid_std_all; plaid_resp_std];
        plaid_cov_all = [plaid_cov_all; plaid_resp_cov];
        
        respPlaid_avg_all = [respPlaid_avg_all; respPlaid_avg];
        respTest_max_all = [respTest_max_all; respTest_max];
        respmask_max_all = [respmask_max_all; respmask_max];

    end

save(fullfile(summaryDir,[svName '_respPlaid_Summary_' area '_' driver '_00' num2str(iexp) '.mat']), 'respPlaid_avg_all', 'respTest_max_all', 'respmask_max_all', 'yfit_all_all_max_rect', 'respplaid_ind_all', 'resp_avg_max_all', 'test_avg_all', 'mask_avg_all', 'yfit_all_all_max', 'plaid_resp_all', 'plaid_std_all', 'plaid_cov_all','p_anova_all_all','b_all_all','amp_all_all','R_square_all_all', 'R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','pha_all_all', 'OSI_all', 'max_dir_all','mouse_list', 'maskPhas_all_all')
%!!!resp plaid script!!!!
    
    
%% Figures comparing areas with SI fit (re-load .mat files generated in first section)
close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC_all','SLC_044','SLC_045','SLC_046','SLC_047', 'SLC_all'); 
area = 'all_areas';
area_list = strvcat('V1','LM','LM', 'LM','LM', 'AL');
narea = length(area_list);
nCells = [];

figure;
    for iA = 1:narea
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];
    
        subplot(2,2,1)
            Rsq_temp = R_square_all_all;
            Rsq_temp(find(Rsq_temp<0)) = 0;
            cdfplot(Rsq_temp(ind,:))
            n = sum(~isnan(Rsq_temp(ind,:)));
            hold on
            xlabel('Rsquared')
            legend(leg_str, 'location', 'southeast')
            
        subplot(2,2,2)
            cdfplot(amp_all_all(ind,:))
            hold on
            xlabel('Phase modulation amplitude')
            legend(leg_str, 'location', 'southeast')
        
        subplot(2,2,3)
            cdfplot(b_all_all(ind))
            hold on
            ylabel('# of cells')
            xlabel('Phase modulation baseline')
            legend(area_list, 'location', 'southeast')          
    end

    print(fullfile(outDir, [svName '_' area '_Compare.pdf']),'-dpdf', '-fillpage') 

    
%% histograms looking at the normalized responses to plaid, test, mask stimuli
close all; clear all; clc;

SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'V1_LM');
svName = 'randPhase';
dateOfAnalysis = '220502';
driver = strvcat('SLC', 'SLC', 'syn'); 
area = 'V1_LM';
area_list = strvcat('V1','LM', 'AL');
narea = length(area_list);

figure;
    for iA = 1:narea

        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary_respPlaid' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        ind = respplaid_ind_all;
 
        leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];
        
%         norm_respPlaid = (amp_all_all(ind)*2) ./ yfit_all_all_max_rect(ind); 
      
        subplot(3,2,1)
%             edges = [0:0.2:1];
            histogram(amp_all_all)
            hold on
             xlim([0 1])
            ylabel('# of cells')
            xlabel('plaid fit amp')
            legend(leg_str)
            
        subplot(3,2,2)
            cdfplot(amp_all_all(ind,:))
            hold on  
            xlabel('phase modulation amplitude')   
    end

        print(fullfile(summaryDir, [svName '_' area '_Compare_respPlaid.pdf']),'-dpdf', '-fillpage')

        

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
dateOfAnalysis = strvcat('220616', '220616', '220829');
driver = strvcat('SLC_a','SLC_a', 'syn_4'); 
area_list = strvcat('V1','LM', 'LM');
narea = length(area_list);
nCells = [];  

figure;
    for iA = 1:narea

        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))    
       
        high_b_ind = find(b_all_all>0.8);
        low_b_ind = find(b_all_all<-0.8);
        
        leg_strH{iA}=[area_list(iA,:) ' n=' num2str(length(high_b_ind))];
        leg_strL{iA}=[area_list(iA,:) ' n=' num2str(length(low_b_ind))];

    SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
    summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
    outDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'V1_LM');
    svName = 'randPhase';
    dateOfAnalysis = '220502';
    driver = strvcat('SLC', 'SLC', 'syn'); 
    area = 'V1_LM';
    area_list = strvcat('V1','LM', 'AL');
    narea = length(area_list);

        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary_respPlaid' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        ind = resp_ind_all;
        
        leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];
        
        subplot(3,2,1)
%         edges = ([0:0.1:1]);
        histogram(amp_all_all(ind))
        hold on
         xlabel('amp')
         ylabel('number of cells')
         ylim([0 175])
         legend(leg_str)
        
        subplot(3,2,2)
%         edges = ([0:.05:0.5]);
        histogram(amp_all_all(high_b_ind))
        hold on
         xlabel('amp')
         ylabel('number of cells')
         title('high b')
         legend(leg_strH)
        
        subplot(3,2,3)
%         edges = ([0:0.05:0.5]);
        histogram(amp_all_all(low_b_ind))
        hold on
         ylim([0 25])
         xlabel('amp')
         ylabel('number of cells')
         title('low b')
         legend(leg_strL)
        
        leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];
        
        subplot(3,2,4)
%         edges = ([-0.1:0.1:2.1]);
        histogram(respPlaid_avg_all(ind))
        hold on
         xlabel('plaid resp-avg-all')
         ylabel('number of cells')
         legend(leg_str)
    end
    
    print(fullfile(outDir, [svName '_' area '_Compare_highBaseline.pdf']),'-dpdf', '-fillpage')
        
%% figure for grant
clear all; close all; clc;

SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'V1_LM');
svName = 'randPhase';
dateOfAnalysis = '220502';
driver = 'SLC'; 
area = 'V1_LM';
area_list = strvcat('V1','LM');
narea = length(area_list);

C = { [.0 .45 .0], [.2 .153 .91]}; %cell array of colors

test_mask_max = [];
    
figure;
    for iA = 1:narea
        
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary' area_list(iA,:) '_' driver '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(3,2,1)
            cdfplot(amp_all_all(ind,:));
            hold on
            xlim([0 .8])
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation amplitude')
            
        subplot(3,2,2)
            cdfplot(b_all_all(ind,:));
            hold on
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation baseline')
            

        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary_respPlaid' area_list(iA,:) '_' driver '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
%      
        test_mask_max = max(test_avg_all, mask_avg_all);
        
        subplot(3,2,3)  
            cdfplot(test_mask_max(ind,:));
            hold on
            legend(leg_str, 'location', 'southeast')
            xlabel('max of test/mask df/f')
        
        subplot(3,2,4)  
            cdfplot(respPlaid_avg_all(ind,:));
            hold on
            legend(leg_str, 'location', 'southeast')
            xlabel('resp plaid avg df/f')

        subplot(3,2,5)  
            cdfplot(resp_avg_max_all(ind,:));
            hold on
            legend(leg_str, 'location', 'southeast')
            xlabel('plaid avg max df/f')
    end
 
    print(fullfile(outDir, [svName '_' area '_GrantFig.pdf']),'-dpdf', '-fillpage')  

    
%% CDF comp of each mouse individually 
clear all; close all; clc;

SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';


figure;

        dateOfAnalysis = strvcat('220614','220614', '220614','220614','220829','220829');
        driver_list = strvcat('SLC_1', 'SLC_2', 'SLC_3', 'SLC_4', 'SLC_5', 'SLC_6'); 
        img_list = strvcat('V1', 'V1', 'V1', 'V1', 'V1', 'V1');
        narea = length(img_list);
        
    for iA = 1:narea
       fprintf([img_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_Summary_' img_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[img_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(2,2,1)
            cdfplot(amp_all_all(ind,:));
            hold on
            xlim([0 .8])
            title('V1 SI fit')
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation amplitude')
            
        subplot(2,2,2)
            cdfplot(b_all_all(ind,:));
            hold on
            title('V1 SI fit')
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation baseline')
        
    end
    
        dateOfAnalysis = strvcat('220623','220623', '220623','220623');
        driver_list = strvcat('SLC_1', 'SLC_2', 'SLC_3', 'SLC_4'); 
        img_list = strvcat('LM', 'LM', 'LM', 'LM');
        narea = length(img_list);
        
    for iA = 1:narea
        fprintf([img_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_Summary_' img_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[img_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(2,2,3)
            cdfplot(amp_all_all(ind,:));
            hold on
            title('LM SI fit')
            xlim([0 .8])
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation amplitude')
            
        subplot(2,2,4)
            cdfplot(b_all_all(ind,:));
            hold on
            title('LM SI fit')
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation baseline')
            end
 movegui('center')

    print(fullfile(outDir, [svName '_SIfit_CompareByAnimal.pdf']),'-dpdf', '-fillpage')  
    
%% CDF comp of each mouse individually for V1
clear all; close all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'V1_LM_V1proj');
svName = 'randPhase';
dateOfAnalysis = strvcat('220614', '220614', '220614', '220614', '220829', '220829', '220829');
driver_list = strvcat('SLC_1', 'SLC_2', 'SLC_3', 'SLC_4', 'SLC_5', 'SLC_6', 'SLC_a'); 
area = 'V1';
img_list = strvcat('V1', 'V1', 'V1', 'V1', 'V1', 'V1', 'V1');
% inj_list = strvcat('V1', 'V1', 'V1', 'V1', 'V1');
narea = length(img_list);

figure;
    for iA = 1:narea
        
        fprintf([img_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_Summary_' img_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[driver_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(2,2,1)
            cdfplot(amp_all_all(ind,:));
            hold on
            xlim([0 .8])
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation amplitude')
            
        subplot(2,2,2)
            cdfplot(b_all_all(ind,:));
            hold on
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation baseline')
                  
    end
 movegui('center')

    print(fullfile(outDir, [svName '_' area '_CompareV1ByAnimal.pdf']),'-dpdf', '-fillpage')  
    
%% CDF comp of each mouse individually for V1 GCaMP8f only
clear all; close all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'V1_LM_V1proj');
svName = 'randPhase';
dateOfAnalysis = strvcat('220616','220616', '220616','220616','220616');
driver_list = strvcat('syn_a', 'syn_1', 'syn_2', 'syn_3', 'syn_4'); 
area = 'V1';
img_list = strvcat('V1', 'V1', 'V1', 'V1', 'V1');
narea = length(img_list);

figure;
    for iA = 1:narea
        
        fprintf([img_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_Summary_' img_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[driver_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(2,2,1)
            cdfplot(amp_all_all(ind,:));
            hold on
            xlim([0 .8])
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation amplitude')
            title('V1')
            
        subplot(2,2,2)
            cdfplot(b_all_all(ind,:));
            hold on
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation baseline')
            title('V1')
                  
    end
 movegui('center')

    print(fullfile(outDir, [svName '_' area '_CompareV1GCaMP8fByAnimal.pdf']),'-dpdf', '-fillpage')  
    
   
%% CDF comp of individual mouse for V1 somas and V1 terminals
clear all; close all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'V1_LM_V1proj');
svName = 'randPhase';
dateOfAnalysis = '220614';
driver_list = strvcat('syn__', 'syn_3'); 
area = 'V1';
img_list = strvcat('V1', 'LM');
inj_list = strvcat('V1', 'V1');
narea = length(img_list);

figure;
    for iA = 1:narea
        fprintf([img_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary_' inj_list(iA,:) '_' img_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[img_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(2,2,1)
            cdfplot(amp_all_all(ind,:));
            hold on
            xlim([0 .8])
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation amplitude')
            
        subplot(2,2,2)
            cdfplot(b_all_all(ind,:));
            hold on
            legend(leg_str, 'location', 'southeast')
            xlabel('Phase modulation baseline')
                  
    end
 movegui('center')

    print(fullfile(outDir, [svName '_' area '_Comparei2722.pdf']),'-dpdf', '-fillpage')  
    
%% compare AL by individual mouse and group avg
clear all; close all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
dateOfAnalysis = strvcat('220623','220623','220829','220829');
driver_list = strvcat('SLC_2', 'SLC_4', 'SLC_8', 'SLC_a'); 
area = 'AL';
area_list = strvcat('AL', 'AL', 'AL', 'AL');
narea = length(area_list);

figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_Summary_' area_list(iA,:) '_' driver_list(iA,:) '.mat'])))
    ind = resp_ind_all;
    leg_str{iA}=[area_list(iA,:) ' ' driver_list(iA,:) ' n=' num2str(length(ind))];
    nCells = length(ind);

    subplot(2,2,1)
        cdfplot(amp_all_all(ind,:));
        hold on
        xlim([0 .8])
        legend(leg_str, 'location', 'southeast')
        xlabel('Phase modulation amplitude')

    subplot(2,2,2)
        cdfplot(b_all_all(ind,:));
        hold on
        legend(leg_str, 'location', 'southeast')
        xlabel('Phase modulation baseline')          
end

%%Pull V1 axon data in AL now (nothing is wrong with the code below, just
%%only wanted somas
% 
% dateOfAnalysis = strvcat('220623','220623','220623','220623','220623');
% driver_list = strvcat('syn_1', 'syn_2', 'syn_3', 'syn_a'); 
% area = 'AL';
% area_list = strvcat('AL', 'AL', 'AL', 'AL');
% narea = length(area_list);
% 
% for iA = 1:narea
%     fprintf([area_list(iA,:) '\n'])
%     load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_Summary_' area_list(iA,:) '_' driver_list(iA,:) '.mat'])))
%     ind = resp_ind_all;
%     leg_str{iA}=[area_list(iA,:) ' ' driver_list(iA,:) ' n=' num2str(length(ind))];
%     nCells = length(ind);
% 
%     subplot(2,2,3)
%         cdfplot(amp_all_all(ind,:));
%         hold on
%         xlim([0 .8])
%         legend(leg_str, 'location', 'southeast')
%         xlabel('Phase modulation amplitude')
% 
%     subplot(2,2,4)
%         cdfplot(b_all_all(ind,:));
%         hold on
%         legend(leg_str, 'location', 'southeast')
%         xlabel('Phase modulation baseline')          
% end
%  movegui('center')
 
    print(fullfile(outDir, [svName '_' area '_Compare.pdf']),'-dpdf', '-fillpage') 

    
%% LM, V1, AL MI yfit amp and baseline with shuffled
close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
dateOfAnalysis = strvcat('220616', '220616','220829');
driver = strvcat('SLC_a','SLC_a', 'SLC_a'); 
area = 'V1_LM_AL';
area_list = strvcat('V1','LM', 'AL');
narea = length(area_list);
nCells = [];

figure;
    for iA = 1:narea
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];
            
        subplot(2,2,1)
            cdfplot(amp_all_all(ind,:))
            hold on
            ylabel('# of cells')
            xlabel('Phase modulation amplitude')
            legend(leg_str, 'location', 'southeast')
        
        subplot(2,2,2)
            cdfplot(b_all_all(ind))
            hold on
            ylabel('# of cells')
            xlabel('Phase modulation baseline')
            legend(area_list, 'location', 'southeast')
    end

    for iA = 1:narea
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];
            
        subplot(2,2,3)
            cdfplot(amp_shuf_all(ind,:))
            hold on
            ylabel('# of cells')
            xlabel('Phase modulation amplitude (SHUFF)')
            legend(leg_str, 'location', 'southeast')
        
        subplot(2,2,4)
            cdfplot(b_shuf_all(ind))
            hold on
            ylabel('# of cells')
            xlabel('Phase modulation baseline (SHUFF)')
            legend(area_list, 'location', 'southeast')
    end    
    
    print(fullfile(outDir, [svName '_' area '_CompareV1LM_Shuff.pdf']),'-dpdf', '-fillpage') 

    


