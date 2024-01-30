%identify cells with OSI > 0.5 on the baseline day
OSI_stat_include = find(mean(norm_dir_resp_stat_concat_rect{pre}(:,1,:),3)>0.05);
OSI_loc_include = find(mean(norm_dir_resp_loc_concat_rect{pre}(:,1,:),3)>0.05);


%%
figure
subplot(2,2,1)
h1=cdfplot(mean(osi_stat{pre}(green_all,:),2,'omitmissing'));
hold on
h2 = cdfplot(mean(osi_stat{post}(green_all,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('Pyr stationary')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI')
ylabel('Cumulative distribution')

subplot(2,2,2)
h1=cdfplot(mean(osi_stat{pre}(red_all,:),2,'omitmissing'));
hold on
h2=cdfplot(mean(osi_stat{post}(red_all,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('SST stationary')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI')
ylabel('Cumulative distribution')

subplot(2,2,3)
h1=cdfplot(mean(osi_loc{pre}(green_all,:),2,'omitmissing'));
hold on
h2=cdfplot(mean(osi_loc{post}(green_all,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('Pyr running')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI')
ylabel('Cumulative distribution')

subplot(2,2,4)
h1=cdfplot(mean(osi_loc{pre}(red_all,:),2,'omitmissing'));
hold on
h2=cdfplot(mean(osi_loc{post}(red_all,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('SST running')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI')
ylabel('Cumulative distribution')
%% delta osi
delta_gOSI_stat = osi_stat{post}-osi_stat{pre};
delta_gOSI_loc = osi_loc{post}-osi_loc{pre};


figure

subplot(2,2,1)
h1=cdfplot(delta_gOSI_stat(green_all));
hold on
h2 = cdfplot(delta_gOSI_loc(green_all));
set(h1, 'Color', 'k');
set(h2, 'Color', 'm');
title('Pyr')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI post-pre')
ylabel('Cumulative distribution')

subplot(2,2,2)
h1=cdfplot(delta_gOSI_stat(red_all));
hold on
h2 = cdfplot(delta_gOSI_loc(red_all));
set(h1, 'Color', 'k');
set(h2, 'Color', 'm');
title('SST')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI post-pre')
ylabel('Cumulative distribution')


subplot(2,2,3)
h1=histogram(delta_gOSI_stat(green_all),20,'FaceColor','k');
hold on
h2 = histogram(delta_gOSI_loc(green_all),20,'FaceColor','m');
xlim([-1.5 1.5])
title('Pyr')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI post-pre')


subplot(2,2,4)
h1=histogram(delta_gOSI_stat(red_all),20,'FaceColor','k');
hold on
h2 = histogram(delta_gOSI_loc(red_all),20,'FaceColor','m');
xlim([-1.5 1.5])
title('SST')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI post-pre')


sgtitle('baseline OSI > 0.5')

%%
figure
subplot(2,2,1)
scatter(osi_stat{pre}(green_ind_concat),delta_gOSI_stat(green_ind_concat))
title('Pyr stat')
xlabel('baseline OSI')
ylabel('delta OSI')

subplot(2,2,2)
scatter(osi_stat{pre}(red_ind_concat),delta_gOSI_stat(red_ind_concat))
title('SST stat')
xlabel('baseline OSI')
ylabel('delta OSI')

subplot(2,2,3)
scatter(osi_loc{pre}(green_ind_concat),delta_gOSI_loc(green_ind_concat))
title('Pyr loc')
xlabel('baseline OSI')
ylabel('delta OSI')

subplot(2,2,4)
scatter(osi_loc{pre}(red_ind_concat),delta_gOSI_loc(red_ind_concat))
title('SST loc')
xlabel('baseline OSI')
ylabel('delta OSI')

%%



dirs_for_plotting=dirs-180;

green_dir_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
red_dir_avrg_stat = cell(1,nd); %same for red
green_dir_se_stat = cell(1,nd); %this will be the se across all green cells
red_dir_se_stat = cell(1,nd); %same for red


    for id = 1:nd
       
        green_dir_avrg_stat{id}=nanmean(nanmean(norm_dir_resp_stat_concat{id}(green_include,:,:),3),1);
        green_std=nanstd(nanmean(norm_dir_resp_stat_concat{id}(green_include,:,:),3),[],1);
        green_dir_se_stat{id}=green_std/sqrt(length(green_include));
        green_dir_avrg_stat{id}=circshift(green_dir_avrg_stat{id},4);
        green_dir_se_stat{id}=circshift(green_dir_se_stat{id},4);
        
        red_dir_avrg_stat{id}=nanmean(nanmean(norm_dir_resp_stat_concat{id}(red_include,:,:),3),1);
        red_std=nanstd(nanmean(norm_dir_resp_stat_concat{id}(red_include,:,:),3),[],1);
        red_dir_se_stat{id}=red_std/sqrt(length(red_include));
        red_dir_avrg_stat{id}=circshift(red_dir_avrg_stat{id},4);
        red_dir_se_stat{id}=circshift(red_dir_se_stat{id},4);
        clear green_std red_std
        
    end
    
    
    
    green_dir_avrg_loc = cell(1,nd); %this will be the average across all green cells - a single line
    red_dir_avrg_loc = cell(1,nd); %same for red
    green_dir_se_loc = cell(1,nd); %this will be the se across all green cells
    red_dir_se_loc = cell(1,nd); %same for red
    
     for id = 1:nd
       
        green_dir_avrg_loc{id}=nanmean(nanmean(norm_dir_resp_loc_concat{id}(green_include,:,:),3),1);
        green_std=nanstd(nanmean(norm_dir_resp_loc_concat{id}(green_include,:,:),3),[],1);
        green_dir_se_loc{id}=green_std/sqrt(length(green_include));
        green_dir_avrg_loc{id}=circshift(green_dir_avrg_loc{id},4);
        green_dir_se_loc{id}=circshift(green_dir_se_loc{id},4);
        
        red_dir_avrg_loc{id}=nanmean(nanmean(norm_dir_resp_loc_concat{id}(red_include,:,:),3),1);
        red_std=nanstd(nanmean(norm_dir_resp_loc_concat{id}(red_include,:,:),3),[],1);
        red_dir_se_loc{id}=red_std/sqrt(length(red_include));
        red_dir_avrg_loc{id}=circshift(red_dir_avrg_loc{id},4);
        red_dir_se_loc{id}=circshift(red_dir_se_loc{id},4);
        clear green_std red_std
        
    end
    
    
    
    figure
    subplot(2,2,1)
    errorbar(dirs_for_plotting,green_dir_avrg_stat{pre},green_dir_se_stat{pre},'k')
    hold on
    errorbar(dirs_for_plotting,green_dir_avrg_stat{post},green_dir_se_stat{post},'b')
    title(['Stationary Pyr'])
    set(gca, 'TickDir', 'out')
    axis square
    box off
    ylabel('dF/F')
    %ylim([-0.01 .2])
    
    subplot(2,2,2)
    errorbar(dirs_for_plotting,red_dir_avrg_stat{pre},red_dir_se_stat{pre},'k')
    hold on
    errorbar(dirs_for_plotting,red_dir_avrg_stat{post},red_dir_se_stat{post},'b')
    title(['Stationary SST'])
    set(gca, 'TickDir', 'out')
    axis square
    box off
    %ylim([-0.01 .2])
    
    subplot(2,2,3)
    errorbar(dirs_for_plotting,green_dir_avrg_loc{pre},green_dir_se_loc{pre},'k')
    hold on
    errorbar(dirs_for_plotting,green_dir_avrg_loc{post},green_dir_se_loc{post},'b')
    title('Running, Pyr')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xlabel('normalized direction')
    ylabel('dF/F')
   % ylim([-0.01 .3])
    
    
    subplot(2,2,4)
    errorbar(dirs_for_plotting,red_dir_avrg_loc{pre},red_dir_se_loc{pre},'k')
    hold on
    errorbar(dirs_for_plotting,red_dir_avrg_loc{post},red_dir_se_loc{post},'b')
    title('Running, SST')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xlabel('normalized direction')
    %ylim([-0.01 .3])
    
    sgtitle(['Normalized direction tuning, rectified, baseline resp > 5% df/f'])
    
%%
mean_NormDir_stat = cell(1,nd);
for id = 1:nd
    mean_NormDir_stat{id}=mean(norm_dir_resp_stat_concat{id},3);
end

elevatedPref=find(mean_NormDir_stat{pre}(:,1)<mean_NormDir_stat{post}(:,1));
suppressedNonPref=find(mean_NormDir_stat{pre}(:,3)>mean_NormDir_stat{post}(:,3));

el_pyr = intersect(elevatedPref, green_ind_concat);
el_sst = intersect(elevatedPref, red_ind_concat);
sup_pyr = intersect(suppressedNonPref, green_ind_concat);
sup_sst = intersect(suppressedNonPref, red_ind_concat);

chancePyr=((length(el_pyr)/length(green_ind_concat))*(length(sup_pyr)/length(green_ind_concat)))*length(green_ind_concat);
chanceSST=((length(el_sst)/length(red_ind_concat))*(length(sup_sst)/length(red_ind_concat)))*length(red_ind_concat);


table([length(el_pyr);length(sup_pyr);length(intersect(el_pyr,sup_pyr));chancePyr], ...
    [length(el_sst);length(sup_sst);length(intersect(el_sst,sup_sst));chanceSST], ...
    'VariableNames',{'Pyr' 'SST'}, 'RowNames',{'Facilitated at pref'  'Suppressed at non pref' 'Intersection' 'Chance intersection'})
%% converting direction to orientation
mean_NormDir_stat = cell(1,nd);
mean_NormDir_loc = cell(1,nd);
for id = 1:nd
    mean_NormDir_stat{id}=mean(norm_dir_resp_stat_concat{id},3);
    mean_NormDir_loc{id}=mean(norm_dir_resp_loc_concat{id},3);
end


normOri_stat= cell(1,nd);
normOri_loc= cell(1,nd);
for id=1:nd
    normOri_stat{id}=mean_NormDir_stat{id}(:,1:4);
    %normOri_stat{id}=((mean_NormDir_stat{id}(:,1:4))+(mean_NormDir_stat{id}(:,5:8)))/2;
    normOri_loc{id}=mean_NormDir_loc{id}(:,1:4);
    %normOri_loc{id}=((mean_NormDir_loc{id}(:,1:4))+(mean_NormDir_loc{id}(:,5:8)))/2;
end
%%

oris_for_plotting=oris-45;

green_dir_avrg_stat = cell(1,nd); %this will be the average across all green cells - a single line
red_dir_avrg_stat = cell(1,nd); %same for red
green_dir_se_stat = cell(1,nd); %this will be the se across all green cells
red_dir_se_stat = cell(1,nd); %same for red


    for id = 1:nd
       
        green_dir_avrg_stat{id}=nanmean(normOri_stat{id}(green_all,:,:),1);
        green_std=nanstd(normOri_stat{id}(green_all,:,:),[],1);
        green_dir_se_stat{id}=green_std/sqrt(length(green_all));
        green_dir_avrg_stat{id}=circshift(green_dir_avrg_stat{id},1);
        green_dir_se_stat{id}=circshift(green_dir_se_stat{id},1);
        
        red_dir_avrg_stat{id}=nanmean(normOri_stat{id}(red_all,:,:),1);
        red_std=nanstd(normOri_stat{id}(red_all,:,:),[],1);
        red_dir_se_stat{id}=red_std/sqrt(length(red_all));
        red_dir_avrg_stat{id}=circshift(red_dir_avrg_stat{id},1);
        red_dir_se_stat{id}=circshift(red_dir_se_stat{id},1);
        clear green_std red_std

         green_dir_avrg_loc{id}=nanmean(normOri_loc{id}(green_all,:,:),1);
        green_std=nanstd(normOri_loc{id}(green_all,:,:),[],1);
        green_dir_se_loc{id}=green_std/sqrt(length(green_all));
        green_dir_avrg_loc{id}=circshift(green_dir_avrg_loc{id},1);
        green_dir_se_loc{id}=circshift(green_dir_se_loc{id},1);
        
        red_dir_avrg_loc{id}=nanmean(normOri_loc{id}(red_all,:,:),1);
        red_std=nanstd(normOri_loc{id}(red_all,:,:),[],1);
        red_dir_se_loc{id}=red_std/sqrt(length(red_all));
        red_dir_avrg_loc{id}=circshift(red_dir_avrg_loc{id},1);
        red_dir_se_loc{id}=circshift(red_dir_se_loc{id},1);
        clear green_std red_std
        
    end
    
    
   
    
    
    figure
    subplot(2,2,1)
    errorbar(oris_for_plotting,green_dir_avrg_stat{pre},green_dir_se_stat{pre},'k')
    hold on
    errorbar(oris_for_plotting,green_dir_avrg_stat{post},green_dir_se_stat{post},'b')
    title(['Stationary Pyr'])
    set(gca, 'TickDir', 'out')
    axis square
    box off
    ylabel('dF/F')
    %ylim([-0.01 .2])
    
    subplot(2,2,2)
    errorbar(oris_for_plotting,red_dir_avrg_stat{pre},red_dir_se_stat{pre},'k')
    hold on
    errorbar(oris_for_plotting,red_dir_avrg_stat{post},red_dir_se_stat{post},'b')
    title(['Stationary SST'])
    set(gca, 'TickDir', 'out')
    axis square
    box off
    %ylim([-0.01 .2])
    
    subplot(2,2,3)
    errorbar(oris_for_plotting,green_dir_avrg_loc{pre},green_dir_se_loc{pre},'k')
    hold on
    errorbar(oris_for_plotting,green_dir_avrg_loc{post},green_dir_se_loc{post},'b')
    title('Running, Pyr')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xlabel('Orientation relative to pref')
    ylabel('dF/F')
   % ylim([-0.01 .3])


    subplot(2,2,4)
    errorbar(oris_for_plotting,red_dir_avrg_loc{pre},red_dir_se_loc{pre},'k')
    hold on
    errorbar(oris_for_plotting,red_dir_avrg_loc{post},red_dir_se_loc{post},'b')
    title('Running, SST')
    set(gca, 'TickDir', 'out')
    axis square
    box off
    xlabel('Orientation relative to pref')
    %ylim([-0.01 .3])
    
    sgtitle(['Normalized orientation tuning'])
    %% get gOSI
gOSI_stat=cell(1,nd);
gOSI_loc=cell(1,nd);
    for id = 1:nd
        gOSI_stat{id}=gOSI(mean_NormDir_stat{id},dirs)';
        gOSI_loc{id}=gOSI(mean_NormDir_loc{id},dirs)';
    end
%% 

green_include = intersect(green_ind_concat,include);
red_include = intersect(red_ind_concat,include);

figure
subplot(2,2,1)
h1=cdfplot(mean(gOSI_stat{pre}(green_include,:),2,'omitmissing'));
hold on
h2 = cdfplot(mean(gOSI_stat{post}(green_include,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('Pyr stationary')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('gOSI')
ylabel('Cumulative distribution')

subplot(2,2,2)
h1=cdfplot(mean(gOSI_stat{pre}(red_include,:),2,'omitmissing'));
hold on
h2=cdfplot(mean(gOSI_stat{post}(red_include,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('SST stationary')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('gOSI')
ylabel('Cumulative distribution')

subplot(2,2,3)
h1=cdfplot(mean(gOSI_loc{pre}(green_include,:),2,'omitmissing'));
hold on
h2=cdfplot(mean(gOSI_loc{post}(green_include,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('Pyr running')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('gOSI')
ylabel('Cumulative distribution')

subplot(2,2,4)
h1=cdfplot(mean(gOSI_loc{pre}(red_include,:),2,'omitmissing'));
hold on
h2=cdfplot(mean(gOSI_loc{post}(red_include,:),2,'omitmissing'));
set(h1, 'Color', 'k');
set(h2, 'Color', 'b');
title('SST running')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('gOSI')
ylabel('Cumulative distribution')

%%
delta_gOSI_stat = gOSI_stat{post}-gOSI_stat{pre};
delta_gOSI_loc = gOSI_loc{post}-gOSI_loc{pre};


figure

subplot(2,2,1)
h1=cdfplot(delta_gOSI_stat(green_include));
hold on
h2 = cdfplot(delta_gOSI_loc(green_include));
set(h1, 'Color', 'k');
set(h2, 'Color', 'm');
title('Pyr')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI post-pre')
ylabel('Cumulative distribution')

subplot(2,2,2)
h1=cdfplot(delta_gOSI_stat(red_include));
hold on
h2 = cdfplot(delta_gOSI_loc(red_include));
set(h1, 'Color', 'k');
set(h2, 'Color', 'm');
title('SST')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI post-pre')
ylabel('Cumulative distribution')


subplot(2,2,3)
h1=histogram(delta_gOSI_stat(green_include),20,'FaceColor','k');
hold on
h2 = histogram(delta_gOSI_loc(green_include),20,'FaceColor','m');
xlim([-1.5 1.5])
title('Pyr')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI post-pre')


subplot(2,2,4)
h1=histogram(delta_gOSI_stat(red_include),20,'FaceColor','k');
hold on
h2 = histogram(delta_gOSI_loc(red_include),20,'FaceColor','m');
xlim([-1.5 1.5])
title('SST')
set(gca, 'TickDir', 'out')
axis square
box off
xlabel('OSI post-pre')


sgtitle('Delta OSI for stationary (black) and running (purple)')