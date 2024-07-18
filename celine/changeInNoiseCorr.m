figure;
subplot(1,2,1)
boxplot([squeeze(R_by_state(1,pre,green_ind_concat)),squeeze(R_by_state(1,post,green_ind_concat))]);
hold on
scatter([1, 2],[squeeze(R_by_state(1,pre,green_ind_concat)),squeeze(R_by_state(1,post,green_ind_concat))],20,'k', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
% xticklabels({'Stationary','Running'})
ylim([-.6 1])
% ylabel(["Mean R value"]) 
% set(gca,'TickDir','out')        
box off

hold on
subplot(1,2,2)
boxplot([squeeze(R_by_state(4,pre,green_ind_concat)),squeeze(R_by_state(4,post,green_ind_concat))]);
hold on
scatter([1, 2],[squeeze(R_by_state(4,pre,green_ind_concat)),squeeze(R_by_state(4,post,green_ind_concat))],20,'k', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)
% xticklabels({'Stationary','Running'})
ylim([-.6 1])
% ylabel(["Mean R value"]) 
% set(gca,'TickDir','out')        
box off

[h1,p1] =ttest(squeeze(R_by_state(1,pre,green_all)),squeeze(R_by_state(1,post,green_all)));
[h2,p2] =ttest(squeeze(R_by_state(4,pre,green_all)),squeeze(R_by_state(4,post,green_all)));
[mean(squeeze(R_by_state(1,pre,green_all))-squeeze(R_by_state(1,post,green_all))),mean(squeeze(R_by_state(4,pre,green_all))-squeeze(R_by_state(4,post,green_all)))]
[p1*2,p2*2]



delta_pre=(squeeze(R_by_state(1,pre,:))-squeeze(R_by_state(4,pre,:)));
delta_post=(squeeze(R_by_state(1,post,:))-squeeze(R_by_state(4,post,:)));

figure
boxplot([delta_pre,delta_post]);
hold on
scatter([1, 2],[delta_pre,delta_post],20,'k', 'MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.25,'jitter', 'on', 'jitterAmount',.1)


[h,p]=ttest(delta_pre, delta_post)