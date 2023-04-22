%% Looking at data quality -- CDFs of average df/f for plaids and STD for all areas
clear all; close all; clc;

SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'timecourses');
svName = 'randPhase';


%% V1, GCaMP6s

dateOfAnalysis = strvcat('221123','221123','221123','221123','221123', '221123','221123');
driver_list = strvcat('SLC_1','SLC_2','SLC_3','SLC_4', 'SLC_5', 'SLC_6', 'SLC_a'); 
area = 'V1';
img_list = strvcat('V1','V1','V1','V1','V1','V1','V1');
inj_list = strvcat('V1','V1','V1','V1','V1','V1','V1');
narea = length(img_list);

figure;
    for iA = 1:narea
        
        fprintf([driver_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_respPlaid_Summary_' img_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area ' ' driver_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(2,2,1)
            cdfplot(plaid_resp_all(ind,:));
            hold on
            title('GCaMP6s plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('Average df/f')
            
        subplot(2,2,2)
            cdfplot(plaid_std_all(ind,:));
            hold on
            title('GCaMP6s plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('standard deviation')      
            
        subplot(2,2,3)
            cdfplot(plaid_cov_all(ind,:));
            hold on
            title('GCaMP6s plaid resp (all)')
            legend(leg_str, 'location', 'northwest')
            xlabel('coefficient of variation, stdev/mean')  
            
        subplot(2,2,4)
            cdfplot(plaid_cov_all(ind,:));
            hold on
            title('GCaMP6s plaid resp (all)')
            xlim([-1 7])
            legend(leg_str, 'location', 'southeast')
            xlabel('coefficient of variation, stdev/mean') 
    end
    
print(fullfile(outDir, ['AveragePlaidResp_V1_gcamp6s.pdf']), '-dpdf','-fillpage')

%% V1, GCaMP8f

dateOfAnalysis = strvcat('221123','221123','221123','221123','221123');
driver_list = strvcat('syn_2','syn_3','syn_4','syn_5','syn_a'); 
area = 'V1';
img_list = strvcat('V1','V1','V1','V1','V1');
inj_list = strvcat('V1','V1','V1','V1','V1');
narea = length(img_list);
    
        for iA = 1:narea
        
        fprintf([driver_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_respPlaid_Summary_' img_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area ' ' driver_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(2,2,1)
            cdfplot(plaid_resp_all(ind,:));
            hold on
            title('GCaMP8f plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('Average df/f')
            
        subplot(2,2,2)
            cdfplot(plaid_std_all(ind,:));
            hold on
            title('GCaPMP8f plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('standard deviation')   
            
        subplot(2,2,3)
            cdfplot(plaid_cov_all(ind,:));
            hold on
            title('GCaMP8f plaid resp (all)')
            legend(leg_str, 'location', 'northwest')
            xlabel('coefficient of variation, stdev/mean')  
            
        subplot(2,2,4)
            cdfplot(plaid_cov_all(ind,:));
            hold on
            title('GCaMP8f plaid resp (all)')
            xlim([-1 7])
            legend(leg_str, 'location', 'southeast')
            xlabel('coefficient of variation, stdev/mean') 
    end
    
 movegui('center')
 print(fullfile(outDir, ['AveragePlaidResp_V1_gcamp8f.pdf']), '-dpdf','-fillpage')
 
 %% LM and AL

dateOfAnalysis = strvcat('221123','221123','221123','221123','221123');
driver_list = strvcat('SLC_1','SLC_2','SLC_3','SLC_4','SLC_a'); 
area = 'LM';
img_list = strvcat('LM','LM','LM','LM','LM');
inj_list = strvcat('LM','LM','LM','LM','LM');
narea = length(img_list);

figure;
    for iA = 1:narea
        
        fprintf([driver_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_respPlaid_Summary_' img_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area ' ' driver_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(3,2,1)
            cdfplot(plaid_resp_all(ind,:));
            hold on
            title('LM plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('Average df/f')
            
        subplot(3,2,2)
            cdfplot(plaid_std_all(ind,:));
            hold on
            title('LM plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('standard deviation')      
    end
    
    
dateOfAnalysis = strvcat('221123','221123','221123','221123');
driver_list = strvcat('SLC_2','SLC_4','SLC_8','SLC_a'); 
area = 'AL';
img_list = strvcat('AL','AL','AL','AL');
inj_list = strvcat('AL','AL','AL','AL');
narea = length(img_list);
    
        for iA = 1:narea
        
        fprintf([driver_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_respPlaid_Summary_' img_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area ' ' driver_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(3,2,3)
            cdfplot(plaid_resp_all(ind,:));
            hold on
            title('AL plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('Average df/f')
            
        subplot(3,2,4)
            cdfplot(plaid_std_all(ind,:));
            hold on
            title('AL plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('standard deviation')     
        
        subplot(3,2,5)
            cdfplot(plaid_cov_all(ind,:));
            hold on
            title('AL plaid resp (all)')
            legend(leg_str, 'location', 'northwest')
            xlabel('coefficient of variation, stdev/mean')  
            
        subplot(3,2,6)
            cdfplot(plaid_cov_all(ind,:));
            hold on
            title('AL plaid resp (all)')
            xlim([-1 7])
            legend(leg_str, 'location', 'southeast')
            xlabel('coefficient of variation, stdev/mean') 
    end
    
 movegui('center')
 print(fullfile(outDir, ['AveragePlaidResp_LM_AL.pdf']), '-dpdf','-fillpage')
 
 %% All areas

dateOfAnalysis = strvcat('221123','221123','221123');
driver_list = strvcat('SLC_a','syn_a','SLC_a'); 
img_list = strvcat('V1','V1','AL');
narea = length(img_list);

figure;
    for iA = 1:narea
        
        fprintf([driver_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([dateOfAnalysis(iA,:) '_' svName '_respPlaid_Summary_' img_list(iA,:) '_' driver_list(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[img_list(iA,:) ' ' driver_list(iA,:) ' n=' num2str(length(ind))];
        nCells = length(ind);
        
        subplot(2,2,1)
            cdfplot(plaid_resp_all(ind,:));
            hold on
            title('plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('Average df/f')
            
        subplot(2,2,2)
            cdfplot(plaid_std_all(ind,:));
            hold on
            title('plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('standard deviation')      
            
        subplot(2,2,3)
            cdfplot(plaid_cov_all(ind,:));
            hold on
            title('plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('CoV (stdev/mean)') 
            
        subplot(2,2,4)
            cdfplot(plaid_cov_all(ind,:));
            hold on
            xlim([-1 7])
            title('plaid resp (all)')
            legend(leg_str, 'location', 'southeast')
            xlabel('CoV (stdev/mean)') 
    end

     print(fullfile(outDir, ['AveragePlaidResp_All.pdf']), '-dpdf','-fillpage')


%% print cell mask



    figure; movegui('center')
    imagesc(mask_cell)
    base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
    outDir = fullfile(base, 'x_temp');
    print(fullfile(outDir, ['cellmask_220620_i2722.pdf']), '-dpdf')
     

    
%% plot scatter of  CoV (std/mean)  by  phase modulation amplitude

close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC_008','SLC_010','SLC_011','SLC_014','SLC_042','SLC_043','SLC_044','SLC_045','SLC_046', 'SLC_007', 'SLC_029', 'SLC_031'); 
area = 'all_exp';
area_list = strvcat('V1','V1', 'V1', 'V1', 'V1', 'V1', 'LM','LM', 'LM', 'AL', 'AL', 'AL');
narea = length(area_list);
nCells = [];

s = 1;
figure;
    for iA = 1:narea
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([svName '_respPlaid_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];
    
        subplot(4,3,s)
            scatter(amp_all_all(ind),plaid_cov_all(ind))
            hold on
            refline(1)
            ylim([0 5])
            xlabel('phase modulation amplitude')
            ylabel('coefficient of variation, stdev/mean')
         s=s+1;      
    end

%     print(fullfile(outDir, [svName '_' area '_CoVbyModulationAmplitude.pdf']),'-dpdf', '-fillpage') 





     