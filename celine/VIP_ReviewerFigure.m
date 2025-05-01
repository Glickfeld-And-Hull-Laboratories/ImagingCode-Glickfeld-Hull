smallSize=find(runningByCondition(:,4,1));
largeSize=find(runningByCondition(:,4,2));
bothSizes=intersect(smallSize,largeSize);
redBothSizes=intersect(bothSizes,red_ind_concat);
%% get mean timecourses

TCmeans=NaN(90,4);
TCsems=NaN(90,4);

TCmeans(:,1)=nanmean(tc_trial_avrg_stat_concat{pre}(:,redBothSizes,4,1),2);
red_std=nanstd(tc_trial_avrg_stat_concat{pre}(:,redBothSizes,4,1),[],2);
TCsems(:,1)=red_std/sqrt(length(redBothSizes));

TCmeans(:,2)=nanmean(tc_trial_avrg_stat_concat{pre}(:,redBothSizes,4,2),2);
red_std=nanstd(tc_trial_avrg_stat_concat{pre}(:,redBothSizes,4,2),[],2);
TCsems(:,2)=red_std/sqrt(length(redBothSizes));

TCmeans(:,3)=nanmean(tc_trial_avrg_loc_concat{pre}(:,redBothSizes,4,1),2);
red_std=nanstd(tc_trial_avrg_loc_concat{pre}(:,redBothSizes,4,1),[],2);
TCsems(:,3)=red_std/sqrt(length(redBothSizes));

TCmeans(:,4)=nanmean(tc_trial_avrg_loc_concat{pre}(:,redBothSizes,4,2),2);
red_std=nanstd(tc_trial_avrg_loc_concat{pre}(:,redBothSizes,4,2),[],2);
TCsems(:,4)=red_std/sqrt(length(redBothSizes));

z=double(nOn)/double(frame_rate);

%create a time axis in seconds
t=1:90;
t=(t-(double(stimStart)-1))/double(frame_rate);
%% plot timecourses
titles = [{'Stationary 20-deg'},{'Stationary fullfield'},{'Running 20-deg'},{'Running fullfield'}]
figure

for iPlot=1:4
    subplot(2,2,iPlot)
        ylim([-.05 .18]);
    hold on
    shadedErrorBar(t,TCmeans(:,iPlot),TCsems(:,iPlot),'k');
    hold on
    line([0,z],[-.01,-.01],'Color','black','LineWidth',2);
    hold on
    line([-1.8,-1.8],[0.01,.06],'Color','black','LineWidth',2);
    title(titles{iPlot})
end

print(fullfile(fnout,['VIPexamples.pdf']),'-dpdf');