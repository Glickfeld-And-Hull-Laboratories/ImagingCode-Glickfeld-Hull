close all; clear all; clc;
doRedChannel = 0;
SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');

doPlot = 1;
ds = ['CrossOriRandPhase_15Hz_ExptList_SG'];
svName = 'randPhase';
eval(ds)
driver = 'SCN';
img_area = {'V1';'L2/3'}; %LM
inj_area = 'V1';

rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);

sigCells = [];
fractSigCells = [];
respCells = [];
resp_ind_all = [];
resptest_ind_all = [];
respmask_ind_all = [];
respplaid_ind_all = [];

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
tmpref_ind=[];
mask_resp_all=[];
test_resp_all=[];

pp_max_all = [];
pp_ind_all = [];

mouse_list = [];

count=0;

totCells = zeros(nexp,1);
start=1;

% V1 L2/3 - 10,11,14,42,43,52,53,54,58,92
% V1 L4 - 93 94

for iexp = [93]
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
        resptest_ind_all = [resptest_ind_all; resptest_ind+sum(totCells(1:iexp-1,:),1)];
        respmask_ind_all = [respmask_ind_all; respmask_ind+sum(totCells(1:iexp-1,:),1)];
        respplaid_ind_all = [respplaid_ind_all; respplaid_ind+sum(totCells(1:iexp-1,:),1)];

        plaid_resp = mean(resp_cell{end,end,1},2);
        mask_resp = mean(resp_cell{end,1,1},2);
        test_resp = mean(resp_cell{1,end,1},2);
        plaid_resp(find(plaid_resp<0)) = 0;
        mask_resp(find(mask_resp<0)) = 0;
        test_resp(find(test_resp<0)) = 0;

        mask_all = resp_cell{end,1,1};
        test_all = resp_cell{1,end,1};

        for ic = 1:nCells
            [h(ic), p(ic)] = ttest2(mask_all(ic,:),test_all(ic,:),'Alpha',0.05,'Tail','both','Vartype','equal');
        end

        tmpref_ind = find(h==1)';
        h=[];
        p=[];

        if start > 1
            tmpref_ind_all = [tmpref_ind_all; count + tmpref_ind];
        else
            tmpref_ind_all = tmpref_ind;
        end
        count=count+nCells;

        mask_resp_all = [mask_resp_all;mask_resp];
        test_resp_all = [test_resp_all;test_resp];

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
        
%         maskPhas_all_all = [maskPhas_all_all; maskPhas_all];
        
        pp_max_all = [pp_max_all; pp_max];
        pp_ind_all = [pp_ind_all; pp_ind];
        
        if exist('resp_avg_max', 'var')
            resp_avg_max_all = [resp_avg_max_all; resp_avg_max];
        else
        end
%     end
start=start+1;
end
save(fullfile(summaryDir,[svName '_Summary_' inj_area '_' driver '_0' num2str(iexp) '.mat']), 'resptest_ind_all','respmask_ind_all', 'respplaid_ind_all', 'pp_max_all', 'pp_ind_all', 'p_anova_all_all','b_all_all','amp_all_all','R_square_all_all', 'R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','pha_all_all', 'OSI_all', 'max_dir_all','mouse_list', 'maskPhas_all_all')
% save(fullfile(summaryDir,[svName '_Summary_' inj_area '_' driver '_all.mat']), 'tmpref_ind_all','mask_resp_all','test_resp_all','resptest_ind_all','respmask_ind_all', 'respplaid_ind_all','pp_max_all', 'pp_ind_all', 'p_anova_all_all','b_all_all','amp_all_all','R_square_all_all', 'R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','pha_all_all', 'OSI_all', 'max_dir_all','mouse_list', 'maskPhas_all_all')


%% summary for resp plaid

close all; clear all; clc;
doRedChannel = 0;
SG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(SG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');

doPlot = 1;
ds = ['CrossOriRandPhase_15Hz_ExptList_SG'];
svName = 'randPhase';
eval(ds)
driver = 'syn';
img_area = {'AL';'L2/3'}; %LM
inj_area = 'AL';

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


for iexp = [74]
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
        
%         maskPhas_all_all = [maskPhas_all_all; maskPhas_all];
        
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
save(fullfile(summaryDir,[svName '_respPlaid_Summary_' inj_area '_' driver '_0' num2str(iexp) '.mat']), 'p_anova_all_all','b_all_all','amp_all_all','R_square_all_all', 'R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','pha_all_all', 'OSI_all', 'max_dir_all','mouse_list', 'maskPhas_all_all')
% save(fullfile(summaryDir,[svName '_respPlaid_Summary_' inj_area '_' driver '_al4.mat']), 'p_anova_all_all','b_all_all','amp_all_all','R_square_all_all', 'R_square_shuf_all','p_anova_shuf_all','b_shuf_all','amp_shuf_all','R_square_shuf_all','yfit_all_all','yfit_shuf_all','plaidSI_all','testPI_all','resp_ind_all','pha_all_all', 'OSI_all', 'max_dir_all','mouse_list', 'maskPhas_all_all', 'mask_avg_all','test_avg_all','respPlaid_avg_all')


%!!!resp plaid script!!!!
    
    
%% Figures comparing areas with SI fit (re-load .mat files generated in first section)
close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC_all','SLC_all','SLC_all'); 
area = 'all_areas';
area_list = strvcat('V1','LM','AL');
narea = length(area_list);
nCells = [];


figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];

    subplot(3,2,1)
        Rsq_temp = R_square_all_all;
        Rsq_temp(find(Rsq_temp<0)) = 0;
        cdfplot(Rsq_temp(ind,:))
        n = sum(~isnan(Rsq_temp(ind,:)));
        hold on
        xlabel('Rsquared')
        legend(leg_str, 'location', 'southeast')

    subplot(3,2,2)
        cdfplot(amp_all_all(ind,:))
        hold on
        xlabel('Phase modulation amplitude')
        legend(leg_str, 'location', 'southeast')

    subplot(3,2,3)
        cdfplot(b_all_all(ind))
        hold on
        ylabel('# of cells')
        xlabel('Phase modulation baseline')
        legend(area_list, 'location', 'southeast')  
        
    subplot(3,2,4)
        cdfplot(amp_shuf_all(ind))
        hold on
        ylabel('# of cells')
        xlabel('Phase modulation amplitude (shuf)')
        xlim([0 0.8])
        legend(area_list, 'location', 'southeast')  
end

stop
    
driver = strvcat('SLC_008','SLC_010','SLC_058','SLC_014','SLC_042','SLC_043','SLC_052','SLC_053','SLC_054','SLC_044','SLC_045','SLC_046','SLC_047','SLC_048','SLC_007','SLC_029','SLC_031','SLC_049','SLC_050','SLC_055','SLC_056','SLC_057');
area_list = strvcat('V1','V1','V1','V1','V1','V1','V1','V1','V1','LM','LM','LM','LM','LM','AL','AL','AL','AL','AL','AL','AL','AL');
narea = length(area_list);

for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];  

    if area_list(iA,:) == 'V1'
        ll = '-';
        c = [0 0.4470 0.7410];
    elseif area_list(iA,:) == 'LM'
        ll = '--';
        c = [0.8500 0.3250 0.0980];
    else 
        ll = ':';
        c = [0.9290 0.6940 0.1250];
    end

    subplot(3,2,5)
        [f,x_cdf] = ecdf(amp_all_all(ind));
        plot(x_cdf,f,'Color',c,'LineWidth',1)
        hold on
        ylabel('# of cells')
        xlabel('Phase modulation amplitude')
%         legend(area_list, 'location', 'southeast')          
end

    % print(fullfile(outDir, [svName '_' area '_Compare.pdf']),'-dpdf', '-fillpage') 

    
%% Comparing AL by TF (with SI fit)
close all; clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';


driver = strvcat('SLC_all','SLC_all','SLC_al4','SLC_al4'); 
area = 'all_areas';
area_list = strvcat('AL', 'V1', 'AL', 'V1');
narea = length(area_list);
nCells = [];

for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];  
    
    if iA == 2
        c = [0.2305 0.6289 0.3281];
    elseif iA == 1
        c = [0.3047 0.6133 0.9453];
    elseif iA == 4
        c = [0.1054 0.4258 0.1836];
    else 
        c = [0.1992 0.4179 0.6563];
    end
    
    subplot(4,4,1)
        [f,x_cdf] = ecdf(R_square_all_all(ind));
        plot(x_cdf,f,'Color',c)
        hold on
        ylabel('CDF')
        xlabel('R squared (MI)')
        legend(leg_str, 'location', 'southeast')

    subplot(4,4,3)
        [f,x_cdf] = ecdf(amp_all_all(ind));
        plot(x_cdf,f,'Color',c)
        hold on
        ylabel('CDF')
        xlabel('Phase modulation amplitude (MI)')

    subplot(4,4,5)
        [f,x_cdf] = ecdf(b_all_all(ind));
        plot(x_cdf,f,'Color',c)
        hold on
        ylabel('CDF')
        xlabel('Phase modulation baseline (MI)')

    subplot(4,4,9)
        [f,x_cdf] = ecdf(pp_max_all(ind));
        plot(x_cdf,f,'Color',c)
        hold on
        ylabel('CDF')
        xlabel('Avg df/f to pref plaid')
end

for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_respPlaid_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))];  
    pref_avg = max(test_avg_all,mask_avg_all); 
    
    if iA == 2
        c = [0.2305 0.6289 0.3281];
    elseif iA == 1 
        c = [0.3047 0.6133 0.9453];
    elseif iA == 4 
        c = [0.1054 0.4258 0.1836];
    else 
        c = [0.1992 0.4179 0.6563];
    end
    
    subplot(4,4,2)
        [f,x_cdf] = ecdf(R_square_all_all(ind));
        plot(x_cdf,f,'Color',c)
        hold on
        ylabel('CDF')
        xlabel('R squared (df/f)')
        legend(leg_str, 'location', 'southeast')

    subplot(4,4,4)
        [f,x_cdf] = ecdf(amp_all_all(ind));
        plot(x_cdf,f,'Color',c)
        hold on
        ylabel('CDF')
        xlabel('Phase modulation amplitude (df/f)')

    subplot(4,4,6)
        [f,x_cdf] = ecdf(b_all_all(ind));
        plot(x_cdf,f,'Color',c)
        hold on
        ylabel('CDF')
        xlabel('Phase modulation baseline (df/f)')

    subplot(4,4,7)
        pref_avg = max(test_avg_all,mask_avg_all);
        [f,x_cdf] = ecdf(pref_avg(ind));
        plot(x_cdf,f,'Color',c)
        hold on
        ylabel('CDF')
        xlabel('Avg df/f to pref test/mask')

    subplot(4,4,8)
        [f,x_cdf] = ecdf(respPlaid_avg_all(ind));
        plot(x_cdf,f,'Color',c)
        hold on
        ylabel('CDF')
        xlabel('Avg df/f across plaids')

    subplot(4,4,11)
        load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        scatter(pref_avg(ind),pp_max_all(ind),1,c)
        hold on
        coeff = polyfit(pref_avg(ind),pp_max_all(ind),1);
        xFit = linspace(min(pref_avg(ind)), max(pref_avg(ind)), 30);
        yFit = polyval(coeff, xFit);
        plot(xFit, yFit, '-','Color',c, 'LineWidth', 1)
        ylabel('Avg df/f to \bfpref\rm plaid')
        xlabel('Avg df/f to pref test/mask')
        title('Grating/plaid preference') 
        
    subplot(4,4,10)
        load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        scatter(pref_avg(ind),respPlaid_avg_all(ind),1,c)
        hold on
        coeff = polyfit(pref_avg(ind),respPlaid_avg_all(ind),1);
        xFit = linspace(min(pref_avg(ind)), max(pref_avg(ind)), 30);
        yFit = polyval(coeff, xFit);
        plot(xFit, yFit, '-','Color',c, 'LineWidth', 1)
        ylabel('Avg df/f to \bfavg\rm plaid')
        xlabel('Avg df/f to pref test/mask')
        title('Grating/plaid preference') 
end
    sgtitle('1 Hz v. 4 Hz')
    print(fullfile(outDir, [svName '_CompareALByTF.pdf']),'-dpdf', '-fillpage')


%% Looking at if invariance correlates with response strength
close all; clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';

driver = strvcat('SLC_all','SLC_all', 'syn_074'); 
area = 'all_areas';
area_list = strvcat('V1', 'AL', 'AL');
narea = length(area_list);
nCells = [];


for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    mod_amp = amp_all_all;  
 
    load(fullfile(summaryDir, ([svName '_respPlaid_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))

    ind = resp_ind_all;
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))]; 
    pref_avg = max(test_avg_all,mask_avg_all);  
    
    if iA == 1
        c = [0.2305 0.6289 0.3281];
    elseif iA == 2
        c = [0.3047 0.6133 0.9453];
    else
        c = [0.7 0.7 0.7];
    end
    
    if iA == 1
    subplot(4,2,1)
        [f,x_cdf] = ecdf(pref_avg(ind));
        plot(x_cdf,f)
        hold on
        [f,x_cdf] = ecdf(respPlaid_avg_all(ind));
        plot(x_cdf,f)
        [f,x_cdf] = ecdf(pp_max_all(ind));
        plot(x_cdf,f)
        ylabel('CDF')
        xlabel('Avg df/f')
        title('V1 - gratings (b), all plaids (r), pref plaid (y)') 
        
    subplot(4,2,3)
        histogram(pp_ind_all,8,'FaceColor',c)
        hold on
        ylabel('num cells')
        xlabel('plaid phase preference')
        title('V1')       
    end
    
    if iA == 2
    subplot(4,2,2)
        [f,x_cdf] = ecdf(pref_avg(ind));
        plot(x_cdf,f)
        hold on
        [f,x_cdf] = ecdf(respPlaid_avg_all(ind));
        plot(x_cdf,f)
        [f,x_cdf] = ecdf(pp_max_all(ind));
        plot(x_cdf,f)
        ylabel('CDF')
        xlabel('Avg df/f')
        title('AL - gratings (b), all plaids (r), pref plaid (y)')  
   
    subplot(4,2,4)
        histogram(pp_ind_all(ind),8,'FaceColor',c)
        hold on
        ylabel('num cells')
        xlabel('plaid phase preference')
        title('AL')    
    end
    
%     subplot(4,2,8)
%         scatter(pref_avg(ind),pp_max_all(ind),5,c)
%         hold on
%         coeff = polyfit(pref_avg(ind),pp_max_all(ind),1);
%         xFit = linspace(min(pref_avg(ind)), max(pref_avg(ind)), 30);
%         yFit = polyval(coeff, xFit);
%         plot(xFit, yFit,'Color',c, 'LineWidth', 1)
%         ylabel('Avg df/f to pref plaid')
%         xlabel('Avg df/f to pref test/mask')
%         title('Grating/plaid preference') 
end

% subplot 5: MI by selectivity index
edges = [0:0.2:1];
for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    
    ind = resp_ind_all;
    MI_amp = amp_all_all(ind);
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))]; 
    load(fullfile(summaryDir, ([svName '_respPlaid_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    test_avg_all(test_avg_all<0) = 0;
    mask_avg_all(mask_avg_all<0) = 0;
    SI_all = (abs((test_avg_all(ind)-mask_avg_all(ind)))./(test_avg_all(ind)+mask_avg_all(ind)));  
        
    if iA == 1
        c = [0.2305 0.6289 0.3281];
    elseif iA == 2
        c = [0.3047 0.6133 0.9453];
    else
        c = [0.7 0.7 0.7];
    end
    
    subplot(4,2,5)
        [N, e, bin] = histcounts(SI_all, edges);
        avgMI = [];
        std_MI = [];
        sem_MI = [];
        avgSI = [];
        std_SI = [];
        sem_SI = [];
        for i = 1:length(edges)
            MI = mean(MI_amp(find(bin==i)));
            stdev = std(MI_amp(find(bin==i)));
            sem = stdev ./ sqrt(length(MI_amp(find(bin==i))));
            avgMI = [avgMI; MI];
            std_MI = [std_MI; stdev];
            sem_MI = [sem_MI; sem];
            SI = mean(SI_all(find(bin==i)));
            stdevS = std(SI_all(find(bin==i)));
            semS = stdevS ./ sqrt(length(SI_all(find(bin==i))));
            avgSI = [avgSI; SI];
            std_SI = [std_SI; stdevS];
            sem_SI = [sem_SI; semS];
        end 
        scatter(avgSI,avgMI,25,c)
        hold on
        text(avgSI(1:(end-1)),avgMI(1:(end-1)),sprintfc('  %d',N))
        errorbar(avgSI,avgMI,sem_MI,'Color',[.7 .7 .7],"LineStyle","none")
        errorbar(avgSI,avgMI,sem_SI,"horizontal",'Color',[.7 .7 .7],"LineStyle","none")
        ylabel('avg MI modulation amp')
        xlabel('SI')
        ylim([0 0.5])
        xlim([0 1])
        xticks([0:0.2:1])
end



% subplot 6: MI by df/f average across all plaids
edges = [0.0001    0.0004    0.0013    0.0046    0.0167    0.0599    0.2154    0.7743    2.7826   10.0000];
for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))

    ind = resp_ind_all;
    MI_amp = amp_all_all(ind);  
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))]; 
    load(fullfile(summaryDir, ([svName '_respPlaid_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    plaid_resp = respPlaid_avg_all(ind);
    plaid_resp(plaid_resp<0) = 0;
    
    if iA == 1
        c = [0.2305 0.6289 0.3281];
    elseif iA == 2
        c = [0.3047 0.6133 0.9453];
    else
        c = [0.7 0.7 0.7];
    end
    
    subplot(4,2,6)
        [N, e, bin] = histcounts(plaid_resp, edges);
        avgMI = [];
        std_MI = [];
        sem_MI = [];
        avgPR = [];
        std_PR = [];
        sem_PR = [];
        for i = 1:length(edges)
            MI = mean(MI_amp(find(bin==i)));
            stdev = std(MI_amp(find(bin==i)));
            sem = stdev ./ sqrt(length(MI_amp(find(bin==i))));
            avgMI = [avgMI; MI];
            std_MI = [std_MI; stdev];
            sem_MI = [sem_MI; sem];
            PR = mean(plaid_resp(find(bin==i)));
            stdevP = std(plaid_resp(find(bin==i)));
            semP = stdevP ./ sqrt(length(plaid_resp(find(bin==i))));
            avgPR = [avgPR; PR];
            std_PR = [std_PR; stdevP];
            sem_PR = [sem_PR; semP];
        end
        scatter(avgPR,avgMI,25,c)
        set(gca,'xscale','log')
        hold on
        text(avgPR(1:(end-1)),avgMI(1:(end-1)),sprintfc('  %d',N))
        errorbar(avgPR,avgMI,sem_MI,'Color',[.7 .7 .7],"LineStyle","none")
        errorbar(avgPR,avgMI,sem_PR,"horizontal",'Color',[.7 .7 .7],"LineStyle","none")
        ylabel('avg MI modulation amp')
        xlabel('avg df/f all plaids')
        xlim([0 10])
        ylim([0 0.5])
end

% subplot 7: MI by df/f average across pref plaid
edges = [0.0001    0.0004    0.0013    0.0046    0.0167    0.0599    0.2154    0.7743    2.7826   10.0000];
for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    
    ind = resp_ind_all;
    MI_amp = amp_all_all(ind);
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))]; 

    prefplaid_resp = pp_max_all(ind);
    prefplaid_resp(prefplaid_resp<0) = 0;
    
    if iA == 1
        c = [0.2305 0.6289 0.3281];
    elseif iA == 2
        c = [0.3047 0.6133 0.9453];
    else
        c = [0.7 0.7 0.7];
    end
    
    subplot(4,2,7)
        [N, e, bin] = histcounts(prefplaid_resp, edges);
        avgMI = [];
        std_MI = [];
        sem_MI = [];
        avgPPR = [];
        std_PPR = [];
        sem_PPR = [];
        for i = 1:length(edges)
            MI = mean(MI_amp(find(bin==i)));
            stdev = std(MI_amp(find(bin==i)));
            sem = stdev ./ sqrt(length(MI_amp(find(bin==i))));
            avgMI = [avgMI; MI];
            std_MI = [std_MI; stdev];
            sem_MI = [sem_MI; sem];
            PPR = mean(prefplaid_resp(find(bin==i)));
            stdevPP = std(prefplaid_resp(find(bin==i)));
            semPP = stdevPP ./ sqrt(length(prefplaid_resp(find(bin==i))));
            avgPPR = [avgPPR; PPR];
            std_PPR = [std_PPR; stdevPP];
            sem_PPR = [sem_PPR; semPP];
        end
        scatter(avgPPR,avgMI,25,c)
        set(gca,'xscale','log')
        hold on
        text(avgPPR(1:(end-1)),avgMI(1:(end-1)),sprintfc('  %d',N))
        errorbar(avgPPR,avgMI,sem_MI,'Color',[.7 .7 .7],"LineStyle","none")
        errorbar(avgPPR,avgMI,sem_PPR,"horizontal",'Color',[.7 .7 .7],"LineStyle","none")
        ylabel('avg MI modulation amp')
        xlabel('avg df/f to pref plaid')
        xlim([0 10])
        ylim([0 0.5])
end


% subplot 8: MI by df/f average across pref grating
edges = [0.0001    0.0004    0.0013    0.0046    0.0167    0.0599    0.2154    0.7743    2.7826   10.0000];
for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    
    ind = resp_ind_all;
    MI_amp = amp_all_all(ind);
    load(fullfile(summaryDir, ([svName '_respPlaid_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))]; 
    prefgrat_resp = max(test_avg_all(ind),mask_avg_all(ind)); 
    prefgrat_resp(prefgrat_resp<0) = 0;
    
    if iA == 1
        c = [0.2305 0.6289 0.3281];
    elseif iA == 2
        c = [0.3047 0.6133 0.9453];
    else
        c = [0.7 0.7 0.7];
    end
    
    subplot(4,2,8)
        [N, e, bin] = histcounts(prefgrat_resp, edges);
        avgMI = [];
        std_MI = [];
        sem_MI = [];
        avgPGR = [];
        std_PGR = [];
        sem_PGR = [];
        for i = 1:length(edges)
            MI = mean(MI_amp(find(bin==i)));
            stdev = std(MI_amp(find(bin==i)));
            sem = stdev ./ sqrt(length(MI_amp(find(bin==i))));
            avgMI = [avgMI; MI];
            std_MI = [std_MI; stdev];
            sem_MI = [sem_MI; sem];
            PGR = mean(prefgrat_resp(find(bin==i)));
            stdevPG = std(prefgrat_resp(find(bin==i)));
            semPG = stdevPG ./ sqrt(length(prefgrat_resp(find(bin==i))));
            avgPGR = [avgPGR; PGR];
            std_PGR = [std_PGR; stdevPG];
            sem_PGR = [sem_PGR; semPG];
        end
        scatter(avgPGR,avgMI,25,c)
        set(gca,'xscale','log')
        hold on
        text(avgPGR(1:(end-1)),avgMI(1:(end-1)),sprintfc('  %d',N))
        errorbar(avgPGR,avgMI,sem_MI,'Color',[.7 .7 .7],"LineStyle","none")
        errorbar(avgPGR,avgMI,sem_PGR,"horizontal",'Color',[.7 .7 .7],"LineStyle","none")
        ylabel('avg MI modulation amp')
        xlabel('avg df/f to pref grating')
        xlim([0 10])
        ylim([0 0.5])
end
    sgtitle('Modulation amplitude variation across responsiveness')
    print(fullfile(outDir, [svName '_CompareInvarianceResponse.pdf']),'-dpdf', '-fillpage')

%% Looking at if baseline correlates with response strength
close all; clear all; clc;
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';


driver = strvcat('SLC_all','SLC_all'); 
area = 'all_areas';
area_list = strvcat('V1', 'AL');
narea = length(area_list);
nCells = [];

% subplot 1: Baseline by selectivity index
edges = [0:0.2:1];
for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;
    base = b_all_all(ind);
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))]; 
    load(fullfile(summaryDir, ([svName '_respPlaid_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    test_avg_all(test_avg_all<0) = 0;
    mask_avg_all(mask_avg_all<0) = 0;
    SI_all = (abs((test_avg_all(ind)-mask_avg_all(ind)))./(test_avg_all(ind)+mask_avg_all(ind)));  
        
    if iA == 1
        c = [0.2305 0.6289 0.3281];
    elseif iA == 2
        c = [0.3047 0.6133 0.9453];
    end
    
    subplot(2,2,1)
        [N, e, bin] = histcounts(SI_all, edges);
        avgB = [];
        std_B = [];
        sem_B = [];
        avgSI = [];
        std_SI = [];
        sem_SI = [];
        for i = 1:length(edges)
            B = mean(base(find(bin==i)));
            stdev = std(base(find(bin==i)));
            sem = stdev ./ sqrt(length(base(find(bin==i))));
            avgB = [avgB; B];
            std_B = [std_B; stdev];
            sem_B = [sem_B; sem];
            SI = mean(SI_all(find(bin==i)));
            stdevS = std(SI_all(find(bin==i)));
            semS = stdevS ./ sqrt(length(SI_all(find(bin==i))));
            avgSI = [avgSI; SI];
            std_SI = [std_SI; stdevS];
            sem_SI = [sem_SI; semS];
        end 
        scatter(avgSI,avgB,50,c,'LineWidth', 2)
        hold on
        text(avgSI(1:(end-1)),avgB(1:(end-1)),sprintfc('  %d',N))
        errorbar(avgSI,avgB,sem_B,'Color',[.7 .7 .7],"LineStyle","none")
        errorbar(avgSI,avgB,sem_SI,"horizontal",'Color',[.7 .7 .7],"LineStyle","none")
        ylabel('avg modulation baseline')
        xlabel('SI')
        xlim([0 1])
        ylim([-1 0.5])
        xticks([0:0.2:1])
end



% subplot 6: MI by df/f average across all plaids
edges = [0.0001    0.0004    0.0013    0.0046    0.0167    0.0599    0.2154    0.7743    2.7826   10.0000];
for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;

    base = b_all_all(ind);  
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))]; 
    load(fullfile(summaryDir, ([svName '_respPlaid_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    plaid_resp = respPlaid_avg_all(ind);
    plaid_resp(plaid_resp<0) = 0;
    
    if iA == 1
        c = [0.2305 0.6289 0.3281];
    elseif iA == 2
        c = [0.3047 0.6133 0.9453];
    end
    
    subplot(2,2,2)
        [N, e, bin] = histcounts(plaid_resp, edges);
        avgB = [];
        std_B = [];
        sem_B = [];
        avgPR = [];
        std_PR = [];
        sem_PR = [];
        for i = 1:length(edges)
            B = mean(base(find(bin==i)));
            stdev = std(base(find(bin==i)));
            sem = stdev ./ sqrt(length(base(find(bin==i))));
            avgB = [avgB; B];
            std_B = [std_B; stdev];
            sem_B = [sem_B; sem];
            PR = mean(plaid_resp(find(bin==i)));
            stdevP = std(plaid_resp(find(bin==i)));
            semP = stdevP ./ sqrt(length(plaid_resp(find(bin==i))));
            avgPR = [avgPR; PR];
            std_PR = [std_PR; stdevP];
            sem_PR = [sem_PR; semP];
        end
        scatter(avgPR,avgB,50,c,'LineWidth', 2)
        set(gca,'xscale','log')
        hold on
        text(avgPR(1:(end-1)),avgB(1:(end-1)),sprintfc('  %d',N))
        errorbar(avgPR,avgB,sem_B,'Color',[.7 .7 .7],"LineStyle","none")
        errorbar(avgPR,avgB,sem_PR,"horizontal",'Color',[.7 .7 .7],"LineStyle","none")
        legend(leg_str)
        ylabel('avg modulation baseline')
        xlabel('avg df/f all plaids')
        xlim([0 10])
        ylim([-1 0.5])
end

% subplot 7: MI by df/f average across pref plaid
edges = [0.0001    0.0004    0.0013    0.0046    0.0167    0.0599    0.2154    0.7743    2.7826   10.0000];
for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;
    
    base = b_all_all(ind);
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))]; 

    prefplaid_resp = pp_max_all(ind);
    prefplaid_resp(prefplaid_resp<0) = 0;
    
    if iA == 1
        c = [0.2305 0.6289 0.3281];
    elseif iA == 2
        c = [0.3047 0.6133 0.9453];
    end
    
    subplot(2,2,3)
        [N, e, bin] = histcounts(prefplaid_resp, edges);
        avgB = [];
        std_B = [];
        sem_B = [];
        avgPPR = [];
        std_PPR = [];
        sem_PPR = [];
        for i = 1:length(edges)
            B = mean(base(find(bin==i)));
            stdev = std(base(find(bin==i)));
            sem = stdev ./ sqrt(length(base(find(bin==i))));
            avgB = [avgB; B];
            std_B = [std_B; stdev];
            sem_B = [sem_B; sem];
            PPR = mean(prefplaid_resp(find(bin==i)));
            stdevPP = std(prefplaid_resp(find(bin==i)));
            semPP = stdevPP ./ sqrt(length(prefplaid_resp(find(bin==i))));
            avgPPR = [avgPPR; PPR];
            std_PPR = [std_PPR; stdevPP];
            sem_PPR = [sem_PPR; semPP];
        end
        scatter(avgPPR,avgB,50,c,'LineWidth', 2)
        set(gca,'xscale','log')
        hold on
        text(avgPPR(1:(end-1)),avgB(1:(end-1)),sprintfc('  %d',N))
        errorbar(avgPPR,avgB,sem_B,'Color',[.7 .7 .7],"LineStyle","none")
        errorbar(avgPPR,avgB,sem_PPR,"horizontal",'Color',[.7 .7 .7],"LineStyle","none")
        ylabel('avg modulation baseline')
        xlabel('avg df/f to pref plaid')
        xlim([0 10])
        ylim([-1 0.5])
end


% subplot 8: MI by df/f average across pref grating
edges = [0.0001    0.0004    0.0013    0.0046    0.0167    0.0599    0.2154    0.7743    2.7826   10.0000];
for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;
    base = b_all_all(ind);
    load(fullfile(summaryDir, ([svName '_respPlaid_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    leg_str{iA}=[area_list(iA,:) ' n=' num2str(length(ind))]; 
    prefgrat_resp = max(test_avg_all(ind),mask_avg_all(ind)); 
    prefgrat_resp(prefgrat_resp<0) = 0;
    
    if iA == 1
        c = [0.2305 0.6289 0.3281];
    elseif iA == 2
        c = [0.3047 0.6133 0.9453];
    end
    
    subplot(2,2,4)
        [N, e, bin] = histcounts(prefgrat_resp, edges);
        avgB = [];
        std_B = [];
        sem_B = [];
        avgPGR = [];
        std_PGR = [];
        sem_PGR = [];
        for i = 1:length(edges)
            B = mean(base(find(bin==i)));
            stdev = std(base(find(bin==i)));
            sem = stdev ./ sqrt(length(base(find(bin==i))));
            avgB = [avgB; B];
            std_B = [std_B; stdev];
            sem_B = [sem_B; sem];
            PGR = mean(prefgrat_resp(find(bin==i)));
            stdevPG = std(prefgrat_resp(find(bin==i)));
            semPG = stdevPG ./ sqrt(length(prefgrat_resp(find(bin==i))));
            avgPGR = [avgPGR; PGR];
            std_PGR = [std_PGR; stdevPG];
            sem_PGR = [sem_PGR; semPG];
        end
        scatter(avgPGR,avgB,50,c,'LineWidth', 2)
        set(gca,'xscale','log')
        hold on
        text(avgPGR(1:(end-1)),avgB(1:(end-1)),sprintfc('  %d',N))
        errorbar(avgPGR,avgB,sem_B,'Color',[.7 .7 .7],"LineStyle","none")
        errorbar(avgPGR,avgB,sem_PGR,"horizontal",'Color',[.7 .7 .7],"LineStyle","none")
        ylabel('avg modulation baseline')
        xlabel('avg df/f to pref grating')
        xlim([0 10])
        ylim([-1 0.5])
end
    sgtitle('Modulation baseline variation across responsiveness')
    print(fullfile(outDir, [svName '_CompareBaselineResponse.pdf']),'-dpdf', '-fillpage')
    
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
driver = strvcat('SLC_all','SLC_all', 'SLC_all'); 
area = 'V1_LM_AL';
area_list = strvcat('V1','LM', 'AL');
narea = length(area_list);
nCells = [];

figure;
    for iA = 1:narea
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];
            
        subplot(3,2,1)
            cdfplot(amp_all_all(ind,:))
            hold on
            ylabel('# of cells')
            xlabel('Phase modulation amplitude')
            set(gca,'TickDir','out'); box off; axis square; grid off
            legend(leg_str, 'location', 'southeast')
        
        subplot(3,2,2)
            cdfplot(b_all_all(ind))
            hold on
            ylabel('# of cells')
            xlabel('Phase modulation baseline')
            legend(area_list, 'location', 'southeast')
            
    end

    for iA = 1:narea
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];
            
        subplot(3,2,5)
            histogram(amp_all_all(ind,:))
            hold on
            ylabel('# of cells')
            xlabel('Phase modulation amplitude')
            set(gca,'TickDir','out'); box off; axis square; grid off
            legend(leg_str, 'location', 'southeast')
    end

    for iA = 1:narea
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];
            
        subplot(3,2,3)
            cdfplot(amp_shuf_all(ind,:))
            hold on
            xlim([0 1])
            ylabel('# of cells')
            xlabel('Phase modulation amplitude (SHUFF)')
            set(gca,'TickDir','out'); box off; axis square; grid off
            legend(leg_str, 'location', 'southeast')
        
        subplot(3,2,4)
            cdfplot(b_shuf_all(ind))
            hold on
            xlim([-1 1])
            ylabel('# of cells')
            xlabel('Phase modulation baseline (SHUFF)')
            legend(area_list, 'location', 'southeast')
    end    
    
    for iA = 1:narea
        fprintf([area_list(iA,:) '\n'])
        load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
        ind = resp_ind_all;
        leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];
            
        subplot(3,2,6)
            histogram(amp_shuf_all(ind,:))
            hold on
            xlim([0 1])
            ylabel('# of cells')
            xlabel('Phase modulation amplitude (SHUFF)')
            set(gca,'TickDir','out'); box off; axis square; grid off
            legend(leg_str, 'location', 'southeast')
        
    end    
    print(fullfile(outDir, [svName '_' area '_CompareV1LM_Shuff.pdf']),'-dpdf', '-fillpage') 

    

%% Looking at if amplitude changes for neurons that prefer test or mask

close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC_all'); 
area = 'all_areas';
area_list = strvcat('V1');
narea = length(area_list);
nCells = [];


figure;
for iA = 1
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];


    mpref_ind=intersect(find(mask_resp_all>test_resp_all),tmpref_ind_all); %Make index of cells where mask response > test response; and difference is significant
    tpref_ind=intersect(find(test_resp_all>mask_resp_all),tmpref_ind_all);

    subplot(2,2,1)
        cdfplot(amp_all_all(ind,:))
        hold on
        cdfplot(amp_all_all(intersect(ind,mpref_ind),:))
        cdfplot(amp_all_all(intersect(ind,tpref_ind),:))
        xlabel('Phase modulation amplitude')
        set(gca,'TickDir','out'); box off; axis square; grid off
        legend(['all, n=' num2str(length(ind))],['mpref, n=' num2str(length(intersect(ind,mpref_ind)))],['tpref, n=' num2str(length(intersect(ind,tpref_ind)))],'location', 'southeast')

end

    print(fullfile(outDir, [svName '_' area '_CompareAmpAcrossTestMaskPrefence.pdf']),'-dpdf', '-fillpage') 


%% Figures comparing layers in V1 with SI fit (re-load .mat files generated in first section)
close all; clear all; clc;

base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\sara';
summaryDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary', 'summaries');
outDir = fullfile(base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');
svName = 'randPhase';
driver = strvcat('SLC_all','SCN_all'); 
area = 'all_areas';
area_list = strvcat('V1','V1');
narea = length(area_list);
nCells = [];


figure;
for iA = 1:narea
    fprintf([area_list(iA,:) '\n'])
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];

    subplot(3,3,1)
        Rsq_temp = R_square_all_all;
        Rsq_temp(find(Rsq_temp<0)) = 0;
        cdfplot(Rsq_temp(ind,:))
        n = sum(~isnan(Rsq_temp(ind,:)));
        hold on
        xlabel('Rsquared')
        set(gca,'TickDir','out'); box off; axis square; grid off
        legend(leg_str, 'location', 'southeast')

    subplot(3,3,2)
        cdfplot(amp_all_all(ind,:))
        hold on
        xlabel('Phase modulation amplitude')
        set(gca,'TickDir','out'); box off; axis square; grid off

    subplot(3,3,3)
        cdfplot(b_all_all(ind))
        hold on
        ylabel('# of cells')
        xlabel('Phase modulation baseline')
        set(gca,'TickDir','out'); box off; axis square; grid off
        
    subplot(3,3,4)
        cdfplot(amp_shuf_all(ind))
        hold on
        ylabel('# of cells')
        xlabel('Phase modulation amplitude (shuf)')
        xlim([0 0.8])
        set(gca,'TickDir','out'); box off; axis square; grid off
end

    
driver = strvcat('SLC_all','SCN_093','SCN_094');
area_list = strvcat('V1','V1','V1');
narea = length(area_list);

for iA = 1:narea
    load(fullfile(summaryDir, ([svName '_Summary_' area_list(iA,:) '_' driver(iA,:) '.mat'])))
    ind = resp_ind_all;
    leg_str{iA}=[area_list(iA,:) ' ' driver(iA,:) ' n=' num2str(length(ind))];  

    if iA == 1
        ll = '-';
        c = [0 0.4470 0.7410];
    elseif iA == 2
        ll = '--';
        c = [0.8500 0.3250 0.0980];
    else 
        ll = ':';
        c = [0.9290 0.6940 0.1250];
    end

    subplot(3,3,5)
        [f,x_cdf] = ecdf(R_square_all_all(ind));
        plot(x_cdf,f,'Color',c,'LineWidth',1)
        hold on
        ylabel('# of cells')
        xlabel('Rsquared')
        set(gca,'TickDir','out'); box off; axis square; grid off
        legend(leg_str, 'location', 'southeast')          
    subplot(3,3,6)
        [f,x_cdf] = ecdf(amp_all_all(ind));
        plot(x_cdf,f,'Color',c,'LineWidth',1)
        hold on
        ylabel('# of cells')
        xlabel('Phase modulation amplitude')
        set(gca,'TickDir','out'); box off; axis square; grid off
    subplot(3,3,7)
        [f,x_cdf] = ecdf(b_all_all(ind));
        plot(x_cdf,f,'Color',c,'LineWidth',1)
        hold on
        ylabel('# of cells')
        xlabel('Phase modulation baseline')
        set(gca,'TickDir','out'); box off; axis square; grid off
     subplot(3,3,8)
        [f,x_cdf] = ecdf(amp_shuf_all(ind));
        plot(x_cdf,f,'Color',c,'LineWidth',1)
        hold on
        ylabel('# of cells')
        xlabel('Phase modulation amp (shuf)')
        set(gca,'TickDir','out'); box off; axis square; grid off
end


    print(fullfile(outDir, [svName '_' area '_CompareLayers.pdf']),'-dpdf', '-fillpage') 
