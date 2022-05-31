close all; clear all; clc;
doRedChannel = 0;
SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');

doPlot = 1;
ds = ['CrossOriRandPhase_15Hz_ExptList'];
svName = 'randPhase';
eval(ds)
dateOfAnalysis = '220502';
driver = 'SLC';
area = 'LM'; %LM


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
for iexp = 7:nexp
    if strcmp(expt(iexp).inj_loc,area) & strcmp(expt(iexp).driver,driver)         
        mouse = expt(iexp).mouse;
        mouse_list = strvcat(mouse_list, mouse);
        date = expt(iexp).date;
        
        if isfield(expt,'copFolder') 
            ImgFolder = expt(iexp).copFolder;
        else
            ImgFolder = expt(iexp).coFolder;
        end
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);
        
        if strcmp(expt(iexp).saveLoc,'sara')
            SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
        elseif strcmp(expt(iexp).saveLoc,'lindsey')
            SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
        end
        
        summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
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
        
        resp_avg_max_all = [resp_avg_max_all; resp_avg_max];

    end
end
 

        

save(fullfile(summaryDir,[dateOfAnalysis '_' svName '_Summary' area '_' driver '.mat']), 'resp_avg_max_all', 'yfit_all_all_max', 'amp_all_all_norm', 'b_all_all_norm', 'plaid_resp', 'p_anova_all_all','b_all_all','amp_all_all','R_square_all_all', 'R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','pha_all_all', 'OSI_all', 'max_dir_all','mouse_list', 'maskPhas_all_all')

%% summary for resp plaid

close all; clear all; clc;
doRedChannel = 0;
SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');

doPlot = 1;
ds = ['CrossOriRandPhase_15Hz_ExptList'];
svName = 'randPhase';
eval(ds)
dateOfAnalysis = '220502';
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
for iexp = 7:nexp
    if strcmp(expt(iexp).inj_loc,area) & strcmp(expt(iexp).driver,driver)         
        mouse = expt(iexp).mouse;
        mouse_list = strvcat(mouse_list, mouse);
        date = expt(iexp).date;
        
        if isfield(expt,'copFolder') 
            ImgFolder = expt(iexp).copFolder;
        else
            ImgFolder = expt(iexp).coFolder;
        end
        nrun = length(ImgFolder);
        run_str = catRunName(cell2mat(ImgFolder), nrun);
        
        if strcmp(expt(iexp).saveLoc,'sara')
            SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
        elseif strcmp(expt(iexp).saveLoc,'lindsey')
            SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
        end
        
        summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
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
        
        respPlaid_avg_all = [respPlaid_avg_all; respPlaid_avg];
        respTest_max_all = [respTest_max_all; respTest_max];
        respmask_max_all = [respmask_max_all; respmask_max];

    end
end
 

        

save(fullfile(summaryDir,[dateOfAnalysis '_' svName '_Summary_respPlaid' area '_' driver '.mat']), 'respPlaid_avg_all', 'respTest_max_all', 'respmask_max_all', 'yfit_all_all_max_rect', 'respplaid_ind_all', 'resp_avg_max_all', 'test_avg_all', 'mask_avg_all', 'yfit_all_all_max', 'plaid_resp', 'p_anova_all_all','b_all_all','amp_all_all','R_square_all_all', 'R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','pha_all_all', 'OSI_all', 'max_dir_all','mouse_list', 'maskPhas_all_all')

    
    
%% Figures comparing areas with SI fit (re-load .mat files generated in first section)
close all; clear all; clc;

SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
dateOfAnalysis = '220502';
driver = 'SLC'; 
area = 'V1_LM';
area_list = strvcat('V1','LM');
narea = length(area_list);
nCells = [];

figure;
    for iA = 1:narea

        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary' area_list(iA,:) '_' driver '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];

    
        subplot(3,2,1)
            Rsq_temp = R_square_all_all;
            Rsq_temp(find(Rsq_temp<0)) = 0;
            cdfplot(Rsq_temp(resp_ind_all,:))
            n = sum(~isnan(Rsq_temp(resp_ind_all,:)));
            hold on
            xlabel('Rsquared')
            xlim([0 1])
            legend(leg_str)
            
        subplot(3,2,2)
            cdfplot(amp_all_all(resp_ind_all,:))
            hold on
            xlabel('Sine Amplitude')
            xlim([0 1])
            legend(leg_str)
        
        subplot(3,2,3)
            edges = [0:0.05:1];
            histogram(amp_all_all(ind),edges)
            hold on
            xlim([0 1])
            ylabel('# of cells')
            xlabel('Phase modulation amplitude')
            legend(area_list)
        
        subplot(3,2,4)
            edges = [-1:0.1:1];
            histogram(b_all_all(ind), edges)
            hold on
            ylabel('# of cells')
            xlabel('Phase modulation baseline')
            legend(area_list)
            
        subplot(3,2,5)
            scatter(b_all_all(ind),amp_all_all(ind), 'o')
            hold on 
            vline(0,'--k')
            xlabel('baseline')
            ylabel('modulation amp')
            legend(area_list)
                  
    end

    print(fullfile(summaryDir, [svName '_' area '_Compare.pdf']),'-dpdf', '-fillpage')  

    
%% histograms looking at the normalized responses to plaid, test, mask stimuli

close all; clear all; clc;

SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
dateOfAnalysis = '220502';
driver = 'SLC'; 
area = 'V1_LM';
area_list = strvcat('V1','LM');
narea = length(area_list);


    
figure;
    for iA = 1:narea

        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary_respPlaid' area_list(iA,:) '_' driver '.mat'])))
        ind = respplaid_ind_all;
 
        leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];
        
        norm_respPlaid = (amp_all_all(ind)*2) ./ yfit_all_all_max_rect(ind); 
      
        subplot(3,2,1)
%             edges = [0:0.2:1];
            histogram(norm_respPlaid)
            hold on
%             xlim([0 1])
            ylabel('# of cells')
            xlabel('plaid fit amp')
            legend(leg_str)
            
        subplot(3,2,2)
            cdfplot(norm_respPlaid(ind,:))
            hold on  
            xlabel('phase modulation amplitude')
            
    end

        print(fullfile(summaryDir, [svName '_' area '_Compare_respPlaid.pdf']),'-dpdf', '-fillpage')


figure;
    for iA = 1:narea

        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary' area_list(iA,:) '_' driver '.mat'])))
    
        high_b_ind = find(b_all_all>0.8);
        low_b_ind = find(b_all_all<-0.8);
        
        leg_strH{iA}=[area_list(iA,:) ' n=' num2str(length(high_b_ind))];
        leg_strL{iA}=[area_list(iA,:) ' n=' num2str(length(low_b_ind))];
        
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary_respPlaid' area_list(iA,:) '_' driver '.mat'])))
        ind = resp_ind_all;
        
        leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];
        
        subplot(3,2,1)
        edges = ([0:0.1:1]);
        histogram(amp_all_all(ind), edges)
        hold on
         xlabel('amp')
         ylabel('number of cells')
         ylim([0 175])
         legend(leg_str)
        
        subplot(3,2,2)
        edges = ([0:.05:0.5]);
        histogram(amp_all_all(high_b_ind),edges)
        hold on
         xlabel('amp')
         ylabel('number of cells')
         title('high b')
         legend(leg_strH)
        

        
        subplot(3,2,3)
        edges = ([0:0.05:0.5]);
        histogram(amp_all_all(low_b_ind),edges)
        hold on
         ylim([0 25])
         xlabel('amp')
         ylabel('number of cells')
         title('low b')
         legend(leg_strL)
         
         
        leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];
        
        subplot(3,2,4)
        edges = ([-0.1:0.1:2.1]);
        histogram(respPlaid_avg_all(ind), edges)
        hold on
         xlabel('plaid resp-avg-all')
         ylabel('number of cells')
         legend(leg_str)
         
%         subplot(3,2,3)
%         scatter(resp_avg_max_all(high_b_ind),mask_avg_all(high_b_ind))
%         hold on
%          xlabel('plaid resp-avg-max')
%          ylabel('mask avg')
%          legend(leg_str,'location','southeast')
%         
%         subplot(3,2,4)
%         scatter(resp_avg_max_all(high_b_ind),test_avg_all(high_b_ind))
%         hold on
%          xlabel('plaid resp-avg-max')
%          ylabel('test avg')
%          legend(leg_str, 'location','southeast')
%          
%         subplot(3,2,5)
%         scatter(respPlaid_avg_all(high_b_ind), mask_avg_all(high_b_ind))
%         hold on
%          xlabel('plaid resp avg')
%          ylabel('mask avg')
%          legend(leg_str)
%          
%         subplot(3,2,6)
%         scatter(respPlaid_avg_all(high_b_ind), test_avg_all(high_b_ind))
%         hold on
%          xlabel('plaid resp avg')
%          ylabel('test avg')
%          legend(leg_str)
        
    end
    
    print(fullfile(summaryDir, [svName '_' area '_Compare_highBaseline.pdf']),'-dpdf', '-fillpage')
        
%         test_ind = [];
%         mask_ind = [];
%         
%         for i = 1:length(test_avg_all)
%             if test_avg_all(i)>-0.01 && test_avg_all(i)<0.01
%             test_ind = ([test_ind; i]);
%             else
%             end
%         end
%         
%         for i = 1:length(mask_avg_all)
%             if mask_avg_all(i)>-0.01 && mask_avg_all(i)<0.01
%             mask_ind = ([mask_ind; i]);
%             else
%             end
%         end
%         
%         int_ind = intersect(test_ind, mask_ind);
%            
%         
%         subplot(4,2,1)
%             edges = [0:0.2:3];
%             histogram(resp_avg_max_all(ind),edges)
%             hold on
%             xlim([0 3])
%             ylabel('# of cells')
%             xlabel('plaid resp-avg-max')
%             legend(leg_str_1)  
% 
%         subplot(4,2,2)
%             edges = [-0.5:0.2:2.75];
%             histogram(yfit_all_all_max(ind),edges)
%             hold on
%             xlim([-0.5 2.75])
%             ylabel('# of cells')
%             xlabel('max amp from yfit')
%             legend(leg_str_1)
%             
%             
%         subplot(4,2,3)
%             edges = [-0.5:0.1:3.5];
%             histogram(test_avg_all(ind),edges)
%             hold on
%             xlim([-0.5 3.5])
%             ylabel('# of cells')
%             xlabel('test avg')
%             legend(area_list)
%             
%         subplot(4,2,4)
%             edges = [-0.5:0.1:2.5];
%             histogram(mask_avg_all(ind),edges)
%             hold on
%             xlim([-0.5 2.5])
%             ylabel('# of cells')
%             xlabel('mask avg')
%             legend(area_list)
%             
%         subplot(4,2,5)
%             scatter(test_avg_all(ind),mask_avg_all(ind), 'o')
%             hold on 
%             xlim([-0.25 3.5])
%             ylim([-0.25 3.5])
%             xlabel('test avg')
%             ylabel('mask avg')
%             legend(area_list)
%        
%         leg_str_2{iA}=[area_list(iA,:) ' n=' num2str(length(test_ind))];
%          
%         subplot(4,2,6)
%             scatter(resp_avg_max_all(test_ind),test_avg_all(test_ind), 'o')
%             hold on 
%             xlim([-0.1 1])
%             xlabel('plaid resp-avg-max')
%             ylabel('-0.01 < test avg < 0.01')            
%             legend(leg_str_2)
% 
%         leg_str_3{iA}=[area_list(iA,:) ' n=' num2str(length(mask_ind))];
%         
%         subplot(4,2,7)
%             scatter(resp_avg_max_all(mask_ind),mask_avg_all(mask_ind), 'o')
%             hold on 
%             xlim([-0.1 1.75])
%             xlabel('plaid resp-avg-max')
%             ylabel('-0.1 < mask avg < 0.01')
%             legend(leg_str_3)
%             
%         leg_str_4{iA}=[area_list(iA,:) ' n=' num2str(length(int_ind))];
%             
%         subplot(4,2,8)
%             scatter(resp_avg_max_all(int_ind),mask_avg_all(int_ind), 'o')
%             hold on 
% %             xlim([0 3.5])
% %             ylim([0 3.5])
%             xlabel('plaid resp-avg-max')
%             ylabel('mask avg of smallest test&mask resp')
%             legend(leg_str_4)

%     print(fullfile(summaryDir, [svName '_' area '_Compare_respAvg.pdf']),'-dpdf', '-fillpage')

    


%% figure for grant
clear all; close all; clc;


SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
dateOfAnalysis = '220502';
driver = 'SLC'; 
area = 'V1_LM';
area_list = strvcat('V1','LM');
narea = length(area_list);

C = { [.0 .45 .0], [.2 .153 .91]}; %cell array of colors


% figure;
%     for iA = 1:narea
% 
%         fprintf([area_list(iA,:) '\n'])
%         load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary' area_list(iA,:) '_' driver '.mat'])))
%         ind = resp_ind_all;
%         leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];
%         nCells = length(ind);
%         
% %         subplot(4,2,1)
% % %             edges = [0:0.04:.8];
% %             histogram(amp_all_all(ind), 'FaceColor', C{iA})
% % %             bincounts = histcounts(amp_all_all(ind),edges);
% % %             bincenters = 0.02:0.04:0.78;
% %             hold on
% %             xlim([0 0.8])
% %             ylabel('# of cells')
% %             xlabel('Phase modulation amplitude')
% %             legend(leg_str)
% % %    
% %         subplot(3,2,2)
% %             plot(bincenters,(bincounts/nCells), 'LineWidth', 2, 'Color', C{iA})
% %             hold on
% %             legend(area_list)
% %             xlabel('Phase modulation amplitude')
% %             ylabel('fraction of cells')
%             
%         subplot(3,2,1)
%             cdfplot(amp_all_all(ind,:));
%             hold on
%             legend(leg_str)
%             xlabel('Phase modulation amplitude')
%         
%     end
 
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
            
%         subplot(3,2,3)
%             cdfplot(amp_all_all(ind,:));
%             hold on
%             xlabel('Phase modulation amplitude')
%             
%         subplot(3,2,4)
%             cdfplot(b_all_all(ind,:));
%             hold on
%             xlabel('Phase modulation baseline')
% 
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis '_' svName '_Summary_respPlaid' area_list(iA,:) '_' driver '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
%         
        test_mask_max = max(test_avg_all, mask_avg_all);
        
%         subplot(3,2,3)
%             cdfplot(test_avg_all(ind,:));
%             hold on
%             legend(leg_str)
%             xlabel('test avg df/f')
%             
%         subplot(3,2,4)
%             cdfplot(mask_avg_all(ind,:));
%             hold on
%             legend(leg_str)
%             xlabel('mask avg df/f')

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
 

    print(fullfile(summaryDir, [svName '_' area '_GrantFig.pdf']),'-dpdf', '-fillpage')  

