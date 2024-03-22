preDataloc=tc_trial_avrg_loc_concat{pre}(stimStart+2:(stimStart+nOn+2),red_all,3);
postDataloc=tc_trial_avrg_loc_concat{post}(stimStart+2:(stimStart+nOn+2),red_all,3);

preDataloc = downsample(preDataloc,8,0);
postDataloc = downsample(postDataloc,8,0);

preDataloc = preDataloc';
postDataloc = postDataloc';

data=horzcat(preDataloc,postDataloc);
cellID=categorical(red_all);
data=table(cellID,data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),data(:,7),data(:,8),'VariableNames',{'cellID','D1B1','D1B2','D1B3','D1B4','D2B1','D2B2','D2B3','D2B4'});

w = table(categorical([1 1 1 1 2 2 2 2 ].'), categorical([1 2 3 4 1 2 3 4].'), 'VariableNames', {'day', 'bin'}); % within-desing

rm = fitrm(data, 'D1B1-D2B4 ~ 1', 'WithinDesign', w);
ranova(rm, 'withinmodel', 'day*bin')