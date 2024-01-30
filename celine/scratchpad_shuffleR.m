%% Fig 4A
%Timecourses during locomotion. We are now selecting the subset of SST
%cells for which both locionary and locomotion bouts coincided with trials
%of that cell's preferred direction at all contrasts on both days.

tc_green_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
tc_red_avrg_loc = cell(1,nd); %same for red
tc_green_se_loc = cell(1,nd); %this will be the se across all green cells
tc_red_se_loc = cell(1,nd); %same for red



for id = 1:nd
    for iCon=1:nCon
        
    tc_green_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,green_all,iCon),2);
    green_std=nanstd(tc_trial_avrg_loc_concat{id}(:,green_all,iCon),[],2);
    tc_green_se_loc{id}(:,iCon)=green_std/sqrt(length(green_all));
    
    tc_red_avrg_loc{id}(:,iCon)=nanmean(tc_trial_avrg_loc_concat{id}(:,red_all,iCon),2);
    red_std=nanstd(tc_trial_avrg_loc_concat{id}(:,red_all,iCon),[],2);
    tc_red_se_loc{id}(:,iCon)=red_std/sqrt(length(red_all));
    
    clear green_std red_std
    end
end

%create a 2-second stimulus marker
z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:(size(tc_green_avrg_loc{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);
positions=[1,2;3,4;5,6];
figure
for iCon = 1:nCon
p1=positions(iCon,1);
p2=positions(iCon,2);

subplot(3,2,p1) 
shadedErrorBar(t,tc_red_avrg_loc{pre}(:,iCon),tc_red_se_loc{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_loc{post}(:,iCon),tc_red_se_loc{post}(:,iCon),'b');
ylim([-.02 .35]);
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2);

if iCon==1
    title(['SST',' n = ', num2str(length(red_all))])
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
end
set(gca,'XColor', 'none','YColor','none')

subplot(3,2,p2) 
ylim([-.02 .35]);
hold on
shadedErrorBar(t,tc_green_avrg_loc{pre}(:,iCon),tc_green_se_loc{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t,tc_green_avrg_loc{post}(:,iCon),tc_green_se_loc{post}(:,iCon),'--b','transparent');
hold on
line([0,z],[-.01,-.01],'Color','black','LineWidth',2)

if iCon==1
    title(['Pyr',' n = ', num2str(length(green_all))])
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
end
set(gca,'XColor', 'none','YColor','none')

x0=5;
y0=0;
width=4;
height=9;
set(gcf,'units','inches','position',[x0,y0,width,height])
set(gca,'XColor', 'none','YColor','none')


end  
print(fullfile(fnout,'Fig_4A_vertical.pdf'),'-dpdf');

%% Figure 4A statistics
%from the full dataframe, extract rows for locionary trials for SST and
%Pyr cells seperately. Use all SST and all Pyr cells.
loc_dfof_summary=dfof_summary(find(dfof_summary.behState=='loc'),:);

SST_loc_dfof = loc_dfof_summary(find(loc_dfof_summary.cellType==1),:);
Pyr_loc_dfof = loc_dfof_summary(find(loc_dfof_summary.cellType==0),:);

%for each cell type, construct a mixed model, with random effects to
%account for having multiple cells per mouse and mutiple measurements per
%cell, and fixed effects for contrast, pre/post day, and their interaction
lme_sst_loc= fitlme(SST_loc_dfof,'dfof~contrast*day+(1|mouseID)+(1|cellID:mouseID)');
lme_pyr_loc= fitlme(Pyr_loc_dfof,'dfof~contrast*day+(1|mouseID)+(1|cellID:mouseID)');


%what is the mean df/f for each contrast on each day, for SST cells
vars = ["dfof"];
factors = ["contrast","day"];
meanScoresByFactor = varfun(@nanmean, ...
                            SST_loc_dfof, ...
                            "InputVariables",vars, ...
                            "GroupingVariables",factors)

%what is the mean df/f for each contrast on each day, for Pyr cells
vars = ["dfof"];
factors = ["contrast","day"];
meanScoresByFactor = varfun(@nanmean, ...
                            Pyr_loc_dfof, ...
                            "InputVariables",vars, ...
                            "GroupingVariables",factors)

% pairwise ttests for dfof response at each contrast for SST cells
[sst_h1, sst_p1]= ttest(pref_responses_loc_concat{pre}(red_ind_concat,1),pref_responses_loc_concat{post}(red_ind_concat,1));
[sst_h2, sst_p2]= ttest(pref_responses_loc_concat{pre}(red_ind_concat,2),pref_responses_loc_concat{post}(red_ind_concat,2));
[sst_h3, sst_p3]= ttest(pref_responses_loc_concat{pre}(red_ind_concat,3),pref_responses_loc_concat{post}(red_ind_concat,3));

%corrected for three tests
sst_pvalues = [(sst_p1*3);(sst_p2*3);(sst_p3*3)];

% pairwise ttests for dfof response at each contrast for Pyr cells
[pyr_h1, pyr_p1]= ttest(pref_responses_loc_concat{pre}(green_ind_concat,1),pref_responses_loc_concat{post}(green_ind_concat,1));
[pyr_h2, pyr_p2]= ttest(pref_responses_loc_concat{pre}(green_ind_concat,2),pref_responses_loc_concat{post}(green_ind_concat,2));
[pyr_h3, pyr_p3]= ttest(pref_responses_loc_concat{pre}(green_ind_concat,3),pref_responses_loc_concat{post}(green_ind_concat,3));

%corrected for three tests
pyr_pvalues = [(pyr_p1*3);(pyr_p2*3);(pyr_p3*3)]
contrasts = cons';

table(contrasts,sst_pvalues,pyr_pvalues)
