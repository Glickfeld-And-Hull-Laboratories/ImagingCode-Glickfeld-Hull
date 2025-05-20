
%% get delta values
for id = 1:nd
    pref_responses_stat_rect{id} = pref_responses_stat_concat{id};
    pref_responses_loc_rect{id} = pref_responses_loc_concat{id};
    pref_responses_stat_rect{id}(pref_responses_stat_rect{id} < 0)=0;
    pref_responses_loc_rect{id}(pref_responses_loc_rect{id} < 0)=0;
end

delta_stat=pref_responses_stat_concat{post}-pref_responses_stat_concat{pre};
norm_delta_stat=(pref_responses_stat_rect{post}-pref_responses_stat_rect{pre})./(pref_responses_stat_rect{post}+pref_responses_stat_rect{pre});
delta_loc=pref_responses_loc_concat{post}-pref_responses_loc_concat{pre};
norm_delta_loc=(pref_responses_loc_rect{post}-pref_responses_loc_rect{pre})./(pref_responses_loc_rect{post}+pref_responses_loc_rect{pre});


%% plot contrast response for stationary cells, not matched for locomotion

conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red

for iDrug = 1:2
    
    thisDrug = find(drug == iDrug-1);
    green_thisDrug = intersect(thisDrug,green_ind_concat);
    red_thisDrug = intersect(thisDrug,red_ind_concat);
        
    conResp_green_avrg_stat{iDrug}=mean(norm_delta_stat(green_thisDrug,:),1,'omitmissing');
    green_std=std(norm_delta_stat(green_thisDrug,:),1,'omitmissing');
    conResp_green_se_stat{iDrug}=green_std/sqrt(length(green_thisDrug));
    
    conResp_red_avrg_stat{iDrug}=mean(norm_delta_stat(red_thisDrug,:),1,'omitmissing');
    red_std=std(norm_delta_stat(red_thisDrug,:),1,'omitmissing');
    conResp_red_se_stat{iDrug}=red_std/sqrt(length(red_thisDrug));
    
    clear green_std red_std
 
end


PEG = 2;
DART=1;

figure
subplot(1,2,1) %for the second day
errorbar(cons,conResp_red_avrg_stat{PEG},conResp_red_se_stat{PEG},'color',"#EA9010");
hold on
errorbar(cons,conResp_red_avrg_stat{DART},conResp_red_se_stat{DART},'b');
title('SST')
% ylabel('Post - pre/post+pre')
% xlabel('contrast') 
xlim([0 1.1])
ylim([-.4 .2])
xticks([.25 .5 1])
set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_stat{PEG},conResp_green_se_stat{PEG},'--','color',"#EA9010");
hold on
errorbar(cons,conResp_green_avrg_stat{DART},conResp_green_se_stat{DART},'--b');
title('Pyr')
% xlabel('contrast') 
% set(gca, 'TickDir', 'out')
xlim([0 1.1])
ylim([-.4 .2])
xticks([.25 .5 1])
box off

x0=5;
y0=5;
width=2.2;
height=1.6;
set(gcf,'units','inches','position',[x0,y0,width,height])

%sgtitle(['Stationary, not fully matched' ])

% green_PEG = intersect(find(drug==0),green_ind_concat);
% green_DART = intersect(find(drug==1),green_ind_concat);
% red_PEG = intersect(find(drug==0),red_ind_concat);
% red_DART = intersect(find(drug==1),red_ind_concat);
% 
% %two-sample t-tests for red, stationary
% [h1, p1]= ttest2(norm_delta_stat(red_PEG,1),norm_delta_stat(red_DART,1));
% [h2, p2]= ttest2(norm_delta_stat(red_PEG,2),norm_delta_stat(red_DART,2));
% [h3, p3]= ttest2(norm_delta_stat(red_PEG,3),norm_delta_stat(red_DART,3));
% %correct for 3 tests
% p1*3
% p2*3
% p3*3
% clear h1 p1 h2 p2 h3 p3
% 
% %two-sample t-tests for green, stationary
% [h1, p1]= ttest2(norm_delta_stat(green_PEG,1),norm_delta_stat(green_DART,1));
% [h2, p2]= ttest2(norm_delta_stat(green_PEG,2),norm_delta_stat(green_DART,2));
% [h3, p3]= ttest2(norm_delta_stat(green_PEG,3),norm_delta_stat(green_DART,3));
% %correct for 3 tests
% p1*3
% p2*3
% p3*3
% clear h1 p1 h2 p2 h3 p3
%% plot contrast response for fully matched cells

conResp_green_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_stat = cell(1,nd); %same for red
conResp_green_se_stat = cell(1,nd); %this will be the se across all green cells
conResp_red_se_stat = cell(1,nd); %same for red

for iDrug = 1:2
    
    thisDrug = find(drug == iDrug-1);
    green_thisDrug = intersect(thisDrug,green_all);
    red_thisDrug = intersect(thisDrug,red_all);
        
    conResp_green_avrg_stat{iDrug}=mean(norm_delta_stat(green_thisDrug,:),1,'omitmissing');
    green_std=std(norm_delta_stat(green_thisDrug,:),1,'omitmissing');
    conResp_green_se_stat{iDrug}=green_std/sqrt(length(green_thisDrug));
    
    conResp_red_avrg_stat{iDrug}=mean(norm_delta_stat(red_thisDrug,:),1,'omitmissing');
    red_std=std(norm_delta_stat(red_thisDrug,:),1,'omitmissing');
    conResp_red_se_stat{iDrug}=red_std/sqrt(length(red_thisDrug));
    
    clear green_std red_std
 
end


PEG = 1;
DART=2;

figure
subplot(1,2,1) %for the second day
errorbar(cons,conResp_red_avrg_stat{PEG},conResp_red_se_stat{PEG},'color',"#EA9010");
hold on
errorbar(cons,conResp_red_avrg_stat{DART},conResp_red_se_stat{DART},'b');
% title('SST')
% ylabel('Post - pre/post+pre')
% xlabel('contrast') 
xlim([0 1.2])
ylim([-.4 .2])
xticks([.25 .5 1])
% set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_stat{PEG},conResp_green_se_stat{PEG},'--','color',"#EA9010");
hold on
errorbar(cons,conResp_green_avrg_stat{DART},conResp_green_se_stat{DART},'--b');
% title('Pyr')
% xlabel('contrast') 
% set(gca, 'TickDir', 'out')
xlim([0 1.2])
ylim([-.4 .2])
xticks([.25 .5 1])
box off
x0=5;
y0=5;
width=2.2;
height=1.6;
set(gcf,'units','inches','position',[x0,y0,width,height])

green_PEG = intersect(find(drug==0),green_all);
green_DART = intersect(find(drug==1),green_all);
red_PEG = intersect(find(drug==0),red_all);
red_DART = intersect(find(drug==1),red_all);

%two-sample t-tests for red, stationary
[h1, p1]= ttest2(norm_delta_stat(red_PEG,1),norm_delta_stat(red_DART,1));
[h2, p2]= ttest2(norm_delta_stat(red_PEG,2),norm_delta_stat(red_DART,2));
[h3, p3]= ttest2(norm_delta_stat(red_PEG,3),norm_delta_stat(red_DART,3));
%correct for 3 tests
p1*3
p2*3
p3*3
clear h1 p1 h2 p2 h3 p3

%two-sample t-tests for green, stationary
[h1, p1]= ttest2(norm_delta_stat(green_PEG,1),norm_delta_stat(green_DART,1));
[h2, p2]= ttest2(norm_delta_stat(green_PEG,2),norm_delta_stat(green_DART,2));
[h3, p3]= ttest2(norm_delta_stat(green_PEG,3),norm_delta_stat(green_DART,3));
%correct for 3 tests
p1*3
p2*3
p3*3
clear h1 p1 h2 p2 h3 p3

% for running
conResp_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
conResp_red_avrg_loc = cell(1,nd); %same for red
conResp_green_se_loc = cell(1,nd); %this will be the se across all green cells
conResp_red_se_loc = cell(1,nd); %same for red

for iDrug = 1:2
    
    thisDrug = find(drug == iDrug-1);
    green_thisDrug = intersect(thisDrug,green_all);
    red_thisDrug = intersect(thisDrug,red_all);
        
    conResp_green_avrg_loc{iDrug}=mean(norm_delta_loc(green_thisDrug,:),1,'omitmissing');
    green_std=std(norm_delta_loc(green_thisDrug,:),1,'omitmissing');
    conResp_green_se_loc{iDrug}=green_std/sqrt(length(green_thisDrug));
    
    conResp_red_avrg_loc{iDrug}=mean(norm_delta_loc(red_thisDrug,:),1,'omitmissing');
    red_std=std(norm_delta_loc(red_thisDrug,:),1,'omitmissing');
    conResp_red_se_loc{iDrug}=red_std/sqrt(length(red_thisDrug));
    
    clear green_std red_std
 
end

figure
subplot(1,2,1) %for the second day
errorbar(cons,conResp_red_avrg_loc{PEG},conResp_red_se_loc{PEG},'color',"#EA9010");
hold on
errorbar(cons,conResp_red_avrg_loc{DART},conResp_red_se_loc{DART},'b');
% title('SST')
% ylabel('Post - pre/post+pre')
% xlabel('contrast') 
xlim([0 1.2])
ylim([-.25 .25])
xticks([.25 .5 1])
% set(gca, 'TickDir', 'out')
box off

subplot(1,2,2) %for the first day
errorbar(cons,conResp_green_avrg_loc{PEG},conResp_green_se_loc{PEG},'--','color',"#EA9010");
hold on
errorbar(cons,conResp_green_avrg_loc{DART},conResp_green_se_loc{DART},'--b');
% title('Pyr')
% xlabel('contrast') 
% set(gca, 'TickDir', 'out')
xlim([0 1.2])
ylim([-.25 .25])
xticks([.25 .5 1])
box off

x0=5;
y0=5;
width=2.2;
height=1.6;
set(gcf,'units','inches','position',[x0,y0,width,height])

%two-sample t-tests for red, running 
[h1, p1]= ttest2(norm_delta_loc(red_PEG,1),norm_delta_loc(red_DART,1));
[h2, p2]= ttest2(norm_delta_loc(red_PEG,2),norm_delta_loc(red_DART,2));
[h3, p3]= ttest2(norm_delta_loc(red_PEG,3),norm_delta_loc(red_DART,3));
%correct for 3 tests
p1*3
p2*3
p3*3
clear h1 p1 h2 p2 h3 p3

%two-sample t-tests for green, running
[h1, p1]= ttest2(norm_delta_loc(green_PEG,1),norm_delta_loc(green_DART,1));
[h2, p2]= ttest2(norm_delta_loc(green_PEG,2),norm_delta_loc(green_DART,2));
[h3, p3]= ttest2(norm_delta_loc(green_PEG,3),norm_delta_loc(green_DART,3));
%correct for 3 tests
p1*3
p2*3
p3*3
clear h1 p1 h2 p2 h3 p3



%% make data table to do ANOVA for norm_delta metric
cellID=(1:nKeep_total)';
drugCol=categorical(drug);
%beh_state_col = repmat('sta',nKeep_total,1);
%highR_col = categorical(logical(noiseCorr_OG_concat{pre}(1,:)>0.5))';
norm_delta_stat_table = array2table(norm_delta_stat,'VariableNames',{'S25','S50','S100'});
labels_table =table(mouseID,cellID,highR_col,drugCol,'VariableNames',{'mouseID' 'cellID' 'highR' 'drugCond'});
norm_delta_summary_stat = [labels_table,norm_delta_stat_table];

%beh_state_col = repmat('loc',nKeep_total,1);
norm_delta_loc_table = array2table(norm_delta_loc,'VariableNames',{'L25','L50','L100'});
labels_table =table(mouseID,cellID,drugCol,'VariableNames',{'mouseID' 'cellID' 'drugCond'});
norm_delta_summary_loc = [labels_table,norm_delta_loc_table];

norm_delta_full =horzcat(norm_delta_summary_stat,norm_delta_loc_table);

norm_delta_full_SST = norm_delta_full(red_all,:);
norm_delta_stat_SST = norm_delta_summary_stat(red_all,:);
norm_delta_stat_pyr = norm_delta_summary_stat(green_all,:);
norm_delta_loc_SST = norm_delta_summary_loc(red_all,:);
norm_delta_loc_Pyr = norm_delta_summary_loc(green_all,:);
norm_delta_stat_unMatchedSST = norm_delta_summary_stat(red_ind_concat,:);
norm_delta_stat_unMatchedPyr = norm_delta_summary_stat(green_ind_concat,:);

norm_deltat_stat_weakCorr=norm_delta_summary_stat(redLow,:);
norm_deltat_stat_strongCorr=norm_delta_summary_stat(redHigh,:);

%% run ANOVA for norm_delta metric
%full model
WithinModel = table(categorical([1 1 1 2 2 2 ]'), categorical([1 2 3 1 2 3]'), 'VariableNames', {'behState', 'contrast'}); % within-design
fullDeltaModel_SST = fitrm(norm_delta_full_SST, 'S25-L100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(fullDeltaModel_SST, 'withinmodel', 'behState*contrast');
% Display the results
disp(ranovatbl);

%stationary model for matched cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
statDeltaModel_SST = fitrm(norm_delta_stat_SST, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(statDeltaModel_SST, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);

%running model for matched cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
locDeltaModel_SST = fitrm(norm_delta_loc_SST, 'L25-L100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(locDeltaModel_SST, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);

%stationary model for un-matched cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
statDeltaModel_SST_unmatched = fitrm(norm_delta_stat_unMatchedSST, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(statDeltaModel_SST_unmatched, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);

%stationary model for un-matched pyr
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
statDeltaModel_PYR_unmatched = fitrm(norm_delta_stat_unMatchedPyr, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(statDeltaModel_PYR_unmatched, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);


%stationary model for matched Pyr cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
statDeltaModel_Pyr = fitrm(norm_delta_stat_pyr, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(statDeltaModel_Pyr, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);


%running model for matched Pyr cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
locDeltaModel_Pyr = fitrm(norm_delta_loc_Pyr, 'L25-L100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(locDeltaModel_Pyr, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);
%% norm_delta contrast value for weakly and strongly correlated cells

red_DART = intersect(red_ind_concat,find(drug));
red_PEG= intersect(red_ind_concat,find(~drug));

highRInds = find(noiseCorr_OG_concat{pre}(1,:)>0.5);
lowRInds = find(noiseCorr_OG_concat{pre}(1,:)<=0.5);
redHigh=intersect(highRInds, red_ind_concat);
redLow=intersect(lowRInds, red_ind_concat);

conResp_redHigh_avrg_stat = cell(1,nd); %this will be the average across all redHigh cells - a single line
conResp_redLow_avrg_stat = cell(1,nd); %same for redLow
conResp_redHigh_se_stat = cell(1,nd); %this will be the se across all redHigh cells
conResp_redLow_se_stat = cell(1,nd); %same for redLow

for iDrug = 1:2
    
    thisDrug = find(drug == iDrug-1);
    lowThisDrug = intersect(thisDrug,redLow);
    highThisDrug = intersect(thisDrug,redHigh);
    length(lowThisDrug)
    length(highThisDrug)

    conResp_redHigh_avrg_stat{iDrug}=mean(norm_delta_stat(highThisDrug,:),1,'omitnan');
    redHigh_std=std(norm_delta_stat(highThisDrug,:),1,'omitnan');
    conResp_redHigh_se_stat{iDrug}=redHigh_std/sqrt(length(highThisDrug));
    
    conResp_redLow_avrg_stat{iDrug}=mean(norm_delta_stat(lowThisDrug,:),1,'omitnan');
    redLow_std=std(norm_delta_stat(lowThisDrug,:),1,'omitnan');
    conResp_redLow_se_stat{iDrug}=redLow_std/sqrt(length(lowThisDrug));
    
    clear redHigh_std redLow_std
 
end

PEG = 1;
DART = 2;
figure
subplot(1,2,1) %
errorbar(cons,conResp_redLow_avrg_stat{PEG},conResp_redLow_se_stat{PEG},'color',"#EA9010");
hold on
errorbar(cons,conResp_redLow_avrg_stat{DART},conResp_redLow_se_stat{DART},'b');
% title(['Weak Corr'])
% xlabel('contrast') 
% ylabel('dF/F, pref ori') 
xlim([0 1.2])
ylim([-.6 .2])
xticks([.25 .5 1])
box off

subplot(1,2,2) %
errorbar(cons,conResp_redHigh_avrg_stat{PEG},conResp_redHigh_se_stat{PEG},'color',"#EA9010");
hold on
errorbar(cons,conResp_redHigh_avrg_stat{DART},conResp_redHigh_se_stat{DART},'b');
% title(['Strong Corr'])
% 
% xlabel('contrast') 
% set(gca, 'TickDir', 'out')
xlim([0 1.2])
ylim([-.6 .2])
xticks([.25 .5 1])
box off

x0=5;
y0=5;
width=2.2;
height=1.6;
set(gcf,'units','inches','position',[x0,y0,width,height])
%% anova for weakly and strongly correlated cells


%stationary model for matched cells
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
corrDeltaModel_SST = fitrm(norm_delta_stat_unMatchedSST, 'S25-S100 ~drugCond*highR', 'WithinDesign', WithinModel);
ranovatbl = ranova(corrDeltaModel_SST, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);

%unpaired t-test for norm_delta vs drug in each corr group at 100% contrast

[h1, p1]= ttest2(norm_delta_stat(intersect(red_PEG,highRInds),3),norm_delta_stat(intersect(red_DART,highRInds),3));


%stationary for weak corr
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
weakCorrDeltaModel_SST = fitrm(norm_deltat_stat_weakCorr, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(weakCorrDeltaModel_SST, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);


%stationary for strong corr
WithinModel = table(categorical([1 2 3]'), 'VariableNames', {'contrast'}); % within-design
strongCorrDeltaModel_SST = fitrm(norm_deltat_stat_strongCorr, 'S25-S100 ~drugCond', 'WithinDesign', WithinModel);
ranovatbl = ranova(strongCorrDeltaModel_SST, 'withinmodel', 'contrast');
% Display the results
disp(ranovatbl);

