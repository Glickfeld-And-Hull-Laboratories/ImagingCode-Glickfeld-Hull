t=1:(size(tc_green_avrg_stat_full{1,1,1},1));
t=(t-(double(stimStart)-1))/double(frame_rate);

t2=1:(size(tc_green_avrg_stat{1,1,1},1));
t2=(t2-(double(stimStart)-1))/double(frame_rate);


figure
subplot(1,2,1) %for the first day

%ylim([-.02 .3]);

shadedErrorBar(t,tc_green_avrg_stat_full{pre}(:,iCon),tc_green_se_stat_full{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_green_avrg_stat_full{post}(:,iCon),tc_green_se_stat_full{post}(:,iCon),'b');
hold on
line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title('HT-')
txt1 = ['n = ', num2str(length(haveRunning_green))];
%txt2 = ['post ',num2str(length(find(~isnan(pref_responses_stat_concat{post}(green_ind_concat)))))];
text(-1.5,-0.03,txt1);
%text(0.75,-0.03,txt2,'Color','b');
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')
hold on
shadedErrorBar(t2,tc_green_avrg_stat{pre}(:,iCon),tc_green_se_stat{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t2,tc_green_avrg_stat{post}(:,iCon),tc_green_se_stat{post}(:,iCon),'--b');
hold on
line([0,.2],[-.02,-.02],'Color','black','LineWidth',2);

subplot(1,2,2) %

shadedErrorBar(t,tc_red_avrg_stat_full{pre}(:,iCon),tc_red_se_stat_full{pre}(:,iCon),'k');
hold on
shadedErrorBar(t,tc_red_avrg_stat_full{post}(:,iCon),tc_red_se_stat_full{post}(:,iCon),'b');
hold on
line([0,2],[-.01,-.01],'Color','black','LineWidth',2);
hold on
line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
title('HT+')
txt1 = ['n = ', num2str(length(haveRunning_red))];
%txt2 = ['post ',num2str(length(find(~isnan(pref_responses_stat_concat{post}(red_ind_concat)))))];
text(-1.5,-0.03,txt1);
%text(0.75,-0.03,txt2,'Color','b');
ylabel('dF/F') 
xlabel('s') 
set(gca,'XColor', 'none','YColor','none')
hold on
shadedErrorBar(t2,tc_red_avrg_stat{pre}(:,iCon),tc_red_se_stat{pre}(:,iCon),'--k');
hold on
shadedErrorBar(t2,tc_red_avrg_stat{post}(:,iCon),tc_red_se_stat{post}(:,iCon),'--b');
hold on
line([0,.2],[-.02,-.02],'Color','black','LineWidth',2);

sgtitle(['stationary, contrast = ' num2str(cons(iCon))])

print(fullfile(fnout,'fullVs200MS_timecourses.pdf'),'-dpdf');
