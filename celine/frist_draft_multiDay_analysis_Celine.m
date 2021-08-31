%%
cd(fn_multi)
fid = fopen('results matched.txt','wt');
%% 1 evoked activity
fprintf(fid,'Looking at evoked activity \n');
%was there a significant difference in visually-evoked responsiveness
% between days 1 and 2?
x = length(intersect(resp_ind_match{1}, resp_ind_match{2}));
d1 = length(resp_ind_match{1});
d2 = length(resp_ind_match{2});
still = length(resp_ind_match{1}) - length(intersect(resp_ind_match{1}, resp_ind_match{2}))
new = length(resp_ind_match{2}) -length(intersect(resp_ind_match{1}, resp_ind_match{2}))

changeTable = table([d1], [d2], [x], [still], [new],'VariableNames',{'Responsive day 1' 'Responsive day 2' 'Resp d1 still on d2' 'Became unresponsive' 'Newly responsive'})
clear x d1 d2 still new
writetable(changeTable,fullfile(fn_multi,'change.csv'),'WriteRowNames',true);


PC_ind_match = find(red_ind_match ==0);

%using the previously determined responsiveness, resp_ind_match
%isolate the responsiveness of PV vs pyr cells
PC_h_match = cell(1,2);
red_h_match = cell(1,2);
for i=1:2
  PC_h_match{i} = h_match{i}(:,:,PC_ind_match);
  red_h_match{i} = h_match{i}(:,:,red_ind_match);
end

PC_resp_match = cell(1,2);
red_resp_match = cell(1,2);

for i=1:2
  PC_resp_match{i} = squeeze(logical(sum(sum(PC_h_match{i},1),2)));
  red_resp_match{i} = squeeze(logical(sum(sum(red_h_match{i},1),2)));
end

%comparing the proportion of cells that were responsive to at least one
%stimulus

fprintf(fid, '\nComparing the proportion of putative pyr cells that were \n responsive to at least one stimulus')
[h p] = ttest(PC_resp_match{1},PC_resp_match{2})

if h ==1
    fprintf(fid, ['\nSignificant, p=' num2str(p)])
else
    fprintf(fid, ['\nNot significant, p=' num2str(p)])
end
fprintf(fid, ['\nMean day 1 ', num2str(mean(PC_resp_match{1}))]);
fprintf(fid, ['\nMean day 2 ', num2str(mean(PC_resp_match{2}))]);


fprintf(fid, '\n \nComparing the proportion PV/SOM cells that were \n responsive to at least one stimulus')
[h p] = ttest(red_resp_match{1},red_resp_match{2})
if h ==1
    fprintf(fid, ['\nSignificant, p=' num2str(p)])
else
    fprintf(fid, ['\nNot significant, p=' num2str(p)])
end
fprintf(fid, ['\nMean day 1 ', num2str(mean(red_resp_match{1}))]);
fprintf(fid, ['\nMean day 2 ', num2str(mean(red_resp_match{2}))]);


% comparing dfof in the response window - I want to determine whether the
% magnitude of the response was different
% not sure that comparing dfof values directly makes sense - maybe I should
% look at the magnitude of change from the baseline condition
%I will look at the average response (averaging over trials of each
%stimulus condition)

%for putative pyr cells and for PV cells, find the mean data
PC_avgResp_match = cell(1,2);
red_avgResp_match = cell(1,2);
for i=1:2
  PC_avgResp_match{i} = squeeze(mean(mean(resp_avg_all{i}(:,:,PC_ind_match,1))));
  red_avgResp_match{i} = squeeze(mean(mean(resp_avg_all{i}(:,:,red_ind_match,1))));
end

fprintf(fid, '\n\nComparing the magnitude (df/f) of the response window for \n putative pyr cells')
[h p] = ttest(PC_avgResp_match{1},PC_avgResp_match{2})
if h ==1
    fprintf(fid, ['\nSignificant, p=' num2str(p)]);
else
    fprintf(fid, ['\nNot significant, p=' num2str(p)]);
end
fprintf(fid, ['\nMean day 1 ', num2str(mean(PC_avgResp_match{1}))]);
fprintf(fid, ['\nMean day 2 ', num2str(mean(PC_avgResp_match{2}))]);



fprintf(fid, '\n\nComparing the magnitude (df/f) of the response window for \n PV/SOM cells')
[h p] = ttest(red_avgResp_match{1},red_avgResp_match{2})
if h ==1
    fprintf(fid, ['\nSignificant, p=' num2str(p)]);
else
    fprintf(fid, ['\nNot significant, p=' num2str(p)]);
end
fprintf(fid, ['\nMean day 1 ', num2str(mean(red_avgResp_match{1}))]);
fprintf(fid, ['\nMean day 2 ', num2str(mean(red_avgResp_match{2}))]);



%% 2 spontaneous activity
%was there a significant different in spontaneous (not visually evoked)
%activity between days 1 and 2
fprintf(fid,'\n\n\nLooking at spontaneous activity \n');

%isolate the non-HT cells
fprintf(fid, '\n\nComparing the fraction of frames with activity for \n putative pyr cells')
PC_match_fractActive = cell(1,2);
for i=1:2
   PC_match_fractActive{i} = fractTimeActive_match{i}(PC_ind_match);
   mean(PC_match_fractActive{i});
end

[h p] = ttest(PC_match_fractActive{1},PC_match_fractActive{2})
if h ==1
    fprintf(fid, ['\nSignificant, p=' num2str(p)]);
else
    fprintf(fid, ['\nNot significant, p=' num2str(p)]);
end
fprintf(fid, ['\nMean day 1 ', num2str(mean(PC_match_fractActive{1}))])
fprintf(fid, ['\nMean day 2 ', num2str(mean(PC_match_fractActive{2}))])


%isolate the red matched cells
red_match_fractActive = cell(1,2);
for i=1:2
   red_match_fractActive{i} = fractTimeActive_match{i}(red_ind_match);
   mean(red_match_fractActive{i});
end

fprintf(fid, '\n\nComparing the fraction of frames with activity for \n PV/SOM cells')
[h p] = ttest(red_match_fractActive{1},red_match_fractActive{2})
if h ==1
    fprintf(fid, ['\nSignificant, p=' num2str(p)]);
else
    fprintf(fid, ['\nNot significant, p=' num2str(p)]);
end
fprintf(fid, ['\nMean day 1 ', num2str(mean(red_match_fractActive{1}))])
fprintf(fid, ['\nMean day 2 ', num2str(mean(red_match_fractActive{2}))])



%comparing spontaneous activity between pyr and PV cells

%% 3 tuning
% was there a significant difference in how many orientations a cells
% responded to? Essentially asking about tuning width. Could look at von
% mises k sharpness value instead
fprintf(fid,'\n\n\nLooking at tuning \n');

%comparing the number of stimuli that each cell was responsive to
PC_nResp_match = cell(1,2);
red_nResp_match = cell(1,2);
for i=1:2
  PC_nResp_match{i} = squeeze(sum(sum(PC_h_match{i},1),2));
  red_nResp_match{i} = squeeze(sum(sum(red_h_match{i},1),2));
end


[h p] = ttest(PC_nResp_match{1},PC_nResp_match{2});
mean(PC_nResp_match{1})
mean(PC_nResp_match{2})

[h p] = ttest(red_nResp_match{1},red_nResp_match{2});
mean(red_nResp_match{1});
mean(red_nResp_match{2});

% the number of stimulus conditions a cell responds to could be influenced
% by the number of orientations or by the number of contrasts, so I want to
% differentiate those

%comparing the number of stimuli that each cell was responsive to
%how many different orientations did each cell respond to, regardless of
%contrast?
PC_nOris_match = cell(1,2);
red_nOris_match = cell(1,2);
for i=1:2
  PC_nOris_match{i} = sum(logical(squeeze(sum(PC_h_match{i},2))),1);
  red_nOris_match{i} = sum(logical(squeeze(sum(red_h_match{i},2))),1);
end


fprintf(fid, '\n\nComparing the number of orientations each putative pyr cell \n responded to, regardless of contrast')
[h p] = ttest(PC_nOris_match{1},PC_nOris_match{2});
if h ==1
    fprintf(fid, ['\nSignificant, p=' num2str(p)]);
else
    fprintf(fid, ['\nNot significant, p=' num2str(p)]);
end
fprintf(fid, ['\nMean day 1 ', num2str(mean(PC_nOris_match{1}))]);
fprintf(fid, ['\nMean day 2 ', num2str(mean(PC_nOris_match{2}))]);




fprintf(fid, '\n\nComparing the number of orientations each PV/SOM cell \n responded to, regardless of contrast')
[h p] = ttest(red_nOris_match{1},red_nOris_match{2});
if h ==1
    fprintf(fid, ['\nSignificant, p=' num2str(p)]);
else
    fprintf(fid, ['\nNot significant, p=' num2str(p)]);
end
fprintf(fid, ['\nMean day 1 ', num2str(mean(red_nOris_match{1}))]);
fprintf(fid, ['\nMean day 2 ', num2str(mean(red_nOris_match{2}))]);



%% Von Mises fit
for day =1:nd
    data_ori_resp{day} = resp_avg_match{day}(:,:,:,1);
    data_ori_resp{day}=max(data_ori_resp{day},0);
end

%loop through each cell
%apply von mises fxn to it

vonMises_output = cell(1,nd);

%The inputs should be the orientations that you presented (theta, 1 x nori) and the average response of each cell to each orientation. 
for day = 1:nd
    vonMises_output{day}=zeros(nCells_match, 6);
    for i = 1:nCells_match
        thisCell = data_ori_resp{day}(:,i);
        thisCell=thisCell';
        [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(oris,thisCell);
        vonMises_output{day}(i,:)=[b_hat, k1_hat, R1_hat,u1_hat,sse,R_square];
    end
end

clear thisCell b_hat k1_hat R1_hat u1_hat sse R_square %get rid of the last iteration of these looping variables
save(fullfile(fn_multi,'vonM_output.mat'),'vonMises_output')

%this is to get the fit itself and find the preferred orientation from that
%fit, if we had provided more orientations
fit_oris = [0:1:179];
y_fits=cell(1,nd);
pref_oris_vonM = cell(1,nd);

for day = 1:nd
    y_fits{day}= zeros(nCells_match, size(fit_oris,2));
    pref_oris_vonM{day} = NaN(nCells_match, 1);
    for i = 1:nCells_match
        in = vonMises_output{day}(i,:);
        b_tmp = in(1);
        k1_tmp = in(2);
        R1_tmp = in(3);
        u1_tmp = in(4);
        if R1_tmp ==0
            y_fits{day}(i,:)= NaN;
            pref_oris_vonM{day}(i)= NaN;
        else
            y_fits{day}(i,:) = b_tmp+R1_tmp.*exp(k1_tmp.*(cos(2.*(deg2rad(fit_oris)-u1_tmp))-1)); 
            pref_oris_vonM{day}(i)= fit_oris(find(y_fits{day}(i,:)==max(y_fits{day}(i,:))));
        end
    end
end
clear b_tmp k1_tmp R1_tmp u1_tmp
save(fullfile(fn_multi,'vonM_fits.mat'),'y_fits')
save(fullfile(fn_multi,'pref_oris_vonM.mat'),'pref_oris_vonM')
prefOriCSV=cell2mat(pref_oris_vonM);
sharpnessCSV=cell2mat(vonMises_output);
sharpnessCSV=sharpnessCSV(:,[2 8]);

writematrix(prefOriCSV,fullfile(fn_multi,'prefOri_vonM.csv'));
writematrix(sharpnessCSV,fullfile(fn_multi,'sharpness.csv'));
%% comparing sharpeness and preferred orientation across days for each cell type
figure;
subplot(2,2,1);
scatter(pref_oris_vonM{1}(PC_ind_match),pref_oris_vonM{2}(PC_ind_match), 'MarkerEdgeColor','g')
hold on;
scatter(pref_oris_vonM{1}(red_ind_match),pref_oris_vonM{2}(red_ind_match), 'MarkerEdgeColor','r')
hold off;
xlabel('D1- pref ori')
ylabel('D2- pref ori')
xlim([0 200])
ylim([0 200])
refline(1)
title(['Orientation preference across days'])

subplot(2,2,2);
scatter(vonMises_output{1}(PC_ind_match,2),vonMises_output{2}(PC_ind_match,2), 'MarkerEdgeColor','g')
hold on;
scatter(vonMises_output{1}(red_ind_match,2),vonMises_output{2}(red_ind_match,2), 'MarkerEdgeColor','r')
hold off;
xlabel('D1- sharpness')
ylabel('D2- sharpness')
%xlim([0 .2])
%ylim([0 .2])
refline(1)
title(['Tuning sharpness across days'])





%%
%how many different contrasts did each cell respond to, regardless of
%orientation? - not sure this is very informative. What I really want to
%know is whether the response was different at each contrast, e.g. stronger
%at weak contrasts
PC_nCons_match = cell(1,2);
red_nCons_match = cell(1,2);
for i=1:2
  PC_nCons_match{i} = sum(logical(squeeze(sum(PC_h_match{i},1))),1);
  red_nCons_match{i} = sum(logical(squeeze(sum(red_h_match{i},1))),1);
end
[h p] = ttest(PC_nCons_match{1},PC_nCons_match{2});
mean(PC_nCons_match{1});
mean(PC_nCons_match{2});

[h p] = ttest(red_nCons_match{1},red_nCons_match{2});
mean(red_nCons_match{1});
mean(red_nCons_match{2});


%% 4 contrast response function
% looking only at cells that were responsive - at their preferred ori (the
% ori that gave the highest response) what did the relationship of dfof to
% contrast look like?
%fprintf(fid,'\n\n\nLooking at contrast sensitivity');

PC_pref_resp_match = cell(1,2);
red_pref_resp_match = cell(1,2);
for i=1:2
  PC_pref_resp_match{i} = resp_avg_match_pref{i}(:,find(red_ind_match==0));
  %PC_pref_resp_match{i} = PC_pref_resp_match{i}(:,all(~isnan(PC_pref_resp_match{i})));   % remove NaN columns
  red_pref_resp_match{i} = resp_avg_match_pref{i}(:,find(red_ind_match==1));
  %red_pref_resp_match{i} = red_pref_resp_match{i}(:,all(~isnan(red_pref_resp_match{i})));   % remove NaN columns
end

%plotting to make sure that I'm looking at the same thing as in the
%multi-panel figure
figure;
errorbar(cons, nanmean(PC_pref_resp_match{1},2), nanstd(PC_pref_resp_match{1},[],2)./sqrt(length(find(red_ind_match==0))),'-o')
hold on
title('match HT-')
ylim([-0.05 0.3])
errorbar(cons, nanmean(PC_pref_resp_match{2},2), nanstd(PC_pref_resp_match{2},[],2)./sqrt(length(find(red_ind_match==0))),'-o')
hold off
%% get the contrast response averaged over PCs or interneurons
PC_mean_pref_resp=cell(1,2);
for day_id = 1:2
    PC_mean_pref_resp{day_id} = nanmean(PC_pref_resp_match{day_id},2);
end

red_mean_pref_resp=cell(1,2);
for day_id = 1:2
    red_mean_pref_resp{day_id} = nanmean(red_pref_resp_match{day_id},2);
end

figure;
errorbar(cons, nanmean(red_pref_resp_match{1},2), nanstd(red_pref_resp_match{1},[],2)./sqrt(length(find(red_ind_match==0))),'-o')
hold on
title('match HT+')
ylim([-0.05 0.3])
errorbar(cons, nanmean(red_pref_resp_match{2},2), nanstd(red_pref_resp_match{2},[],2)./sqrt(length(find(red_ind_match==0))),'-o')
hold off
%% get the c50 contrast sensitivity for preferred ori of responsive cells

%% get the c50 contrast sensitivity for each cell, averaging over orientations
con_resp_curves_all = cell(1,2);
c50_all = cell(1,2);
for day = 1:2
    nCells = nC_all(day);
    for iCell = 1:nCells
        iCurve = mean(resp_avg_all{day}(:,:,iCell,1),1);
        con_resp_curves_all{day}(:,iCell)=iCurve;
        c50_all{day}(:,iCell)=find_c50(iCurve,cons);
    end
end

%enter line to save the c50 output

figure;
subplot(1,2,1)
hist(c50_all{1});
title('Distribution fo c50 values Day 1');
vline(nanmean(c50_all{1}));
subplot(1,2,2)
hist(c50_all{2});
title('Distribution of c50 values Day 2');
vline(nanmean(c50_all{2}));

figure;
subplot(1,2,1)
plot(cons,con_resp_curves_all{1},'LineWidth',.005,'color',[.2 .2 .2]);
ylim([-.04 .05])
hold on
plot(cons,mean(con_resp_curves_all{1},2), 'LineWidth',2, 'color','b');
hold on 
vline(nanmean(c50_all{1}));
title('Contrast responses 30-degree');
subplot(1,2,2)
plot(cons,con_resp_curves_all{2},'LineWidth',.005,'color',[.2 .2 .2]);
ylim([-.04 .05])
hold on
plot(cons,mean(con_resp_curves_all{2},2), 'LineWidth',2, 'color','b');
hold on 
vline(nanmean(c50_all{2}));
title('Contrast responses fullfield');
sgtitle('All cells, averaged over orientations');


%% select c50 for matched cells only 
con_resp_curves_match = cell(1,2);
c50_match = cell(1,2);
for day = 1:2
    nCells = nC_match;
    for iCell = 1:nCells
        for iOri = 1:nOri
            iCurve = resp_avg_match{day}(iOri,:,iCell,1);
            con_resp_curves_match{day}(iOri,:,iCell)=iCurve;
            c50_match{day}(iOri,iCell)=find_c50(iCurve,cons);
        end
    end
end
save(fullfile(fnout,'c50_match'),'c50_match');
save(fullfile(fnout,'con_resp_curves_match'),'con_resp_curves_match');

figure;
scatter(nanmean(c50_match{1},1),nanmean(c50_match{2},1))
hold on
scatter(nanmean(c50_match{1}(:,red_ind_match),1),nanmean(c50_match{2}(:,red_ind_match),1), 'filled','r')
xlabel('30Deg C50');
ylabel('FF C50');
xlim([0 1]);
ylim([0 1]);
refline(1);
hold off


c50_match_table = cell2mat(c50_match)';
c50_match_table = reshape(c50_match_table,[(70*8),1]);
%change subsequent parts to reflect reshaping
%move the data to column 2
c50_match_table(:,2)=c50_match_table(:,1);
%enter cell ID as column 1
c50_match_table(:,1)=repmat(1:nC_match,1,16)';
%enter day ID as column 3
c50_match_table(:,3)= repelem(repmat(1:2,1,8),nC_match)';
%enter orientation as column 4
c50_match_table(:,4)= repelem(oris,70)';

%enter cell type as column 5
for i =1:size(c50_match_table,1)
 c50_match_table(i,5)=red_ind_match(c50_match_table(i,1));
end

%make it into a table
c50_match_table = array2table(c50_match_table);
c50_match_table.Properties.VariableNames(1:5) = {'cellID','c50','day','ori','cell_type'};
c50_match_table.cell_type = categorical( c50_match_table.cell_type,0:1,{'PYR' 'IN'});
writetable(c50_match_table,fullfile(fn_multi,'c50_match_table.csv'));


figure;
subplot(1,2,1)
hist(nanmean(c50_match{1},1));
title('Distribution of c50 values 30-degree');
vline(nanmean(nanmean(c50_match{1},1)));
xlim([0.2 .8]);
subplot(1,2,2)
hist(nanmean(c50_match{2},1));
title('Distribution of c50 values fullfield');
vline(nanmean(nanmean(c50_match{2},1)));
xlim([0.2 .8]);
%sgtitle('Matched cells, averaged over orientations');

% Need to update to reflect orientation dimension in c50 data
%t-test for signigicant difference in c50 values
[h p] = ttest(c50_match{1},c50_match{2})
nanmean(c50_match{1})
nanmean(c50_match{2})

%looking at change in c50 for matched cells
change_c50 = NaN(1,nCells_match);
for iCell = 1:nCells_match
   change_c50(iCell) = c50_match{1}(:,iCell)-c50_match{2}(:,iCell);
end

figure;
hist(change_c50)
[h p] = ttest(change_c50)


figure;
subplot(1,2,1)
plot(cons,squeeze(nanmean(con_resp_curves_match{1},1)),'LineWidth',.005,'color',[.2 .2 .2]);
ylim([-.04 .05])
xlim([0.1 1]);
hold on
plot(cons,mean(squeeze(nanmean(con_resp_curves_match{1},1)),2), 'LineWidth',2, 'color','b');
%hold on 
%vline(nanmean(nanmean(c50_match{1},1)));
title('Contrast responses 30-degree');
subplot(1,2,2)
plot(cons,squeeze(nanmean(con_resp_curves_match{2},1)),'LineWidth',.005,'color',[.2 .2 .2]);
ylim([-.04 .05])
xlim([0.1 1]);
hold on
plot(cons,mean(squeeze(nanmean(con_resp_curves_match{2},1)),2), 'LineWidth',2, 'color','b');
%hold on 
%vline(nanmean(nanmean(c50_match{2},1)));
title('Contrast responses fullfield');
%sgtitle('Matched cells, averaged over orientations');


%% convert contrast response to slopes
PC_match_slopes =cell(1,2);
red_match_slopes =cell(1,2);
for i = 1:2
    for iCell = 1:size(PC_pref_resp_match{i},2)
        PC_match_slopes{i}(iCell)=(PC_pref_resp_match{i}(5,iCell)-PC_pref_resp_match{i}(1,iCell))/5;
    end
    for iCell = 1:size(red_pref_resp_match{i},2)
        red_match_slopes{i}(iCell)=(red_pref_resp_match{i}(5,iCell)-red_pref_resp_match{i}(1,iCell))/5;
    end
end
% do an anova/regression for multiple levels of contrast


%fprintf(fid,'\n\nComparing slope of df/f vs. contrast at preferred ori\n for putative pyr cells that were responsive');
%[h p] = ttest(PC_match_slopes{1}, PC_match_slopes{2});

% if h ==1
%     fprintf(fid, ['\nSignificant, p=' num2str(p)]);
% else
%     fprintf(fid, ['\nNot significant, p=' num2str(p)]);
% end
% fprintf(fid, ['\nMean day 1 ', num2str(nanmean(PC_match_slopes{1}))]);
% fprintf(fid, ['\nMean day 1 ', num2str(nanmean(PC_match_slopes{2}))]);
% 
% fprintf(fid,'\n\nComparing slop of df/f vs. contrast at preferred ori\n for PV/SOM cells that were responsive');
% [h, p] = ttest(red_match_slopes{1}, red_match_slopes{2});
% if h ==1
%     fprintf(fid, ['\nSignificant, p=' num2str(p)]);
% else
%     fprintf(fid, ['\nNot significant, p=' num2str(p)]);
% end
% fprintf(fid, ['\nMean day 1 ', num2str(nanmean(red_match_slopes{1}))]);
% fprintf(fid, ['\nMean day 1 ', num2str(nanmean(red_match_slopes{2}))]);
%% looking at potential change in preferred orientation from day 1 to day 2
prefOri_matchb = prefOri_match;
for x = 1:size(prefOri_match{1},2)
    if prefOri_match{1}(x)==0
       prefOri_match{1}(x) = NaN 
    else
       prefOri_match{1}(x) = oris(prefOri_match{1}(x))
    end
end
for x = 1:size(prefOri_match{2},2)
    if prefOri_match{2}(x)==0
       prefOri_match{2}(x) = NaN 
    else
       prefOri_match{2}(x) = oris(prefOri_match{2}(x))
    end
end
figure
subplot(1,2,1)
histogram(prefOri_match{1}(~isnan(prefOri_match{1}(:))),8)
subplot(1,2,2)
histogram(prefOri_match{2}(~isnan(prefOri_match{2}(:))),8)
sgtitle('distribution of preferred orientations day 1 and say 2');

prefOri_diff = NaN(1,size(prefOri_match{1},2));
for x = 1:size(prefOri_match{1},2)
    prefOri_diff(x) = prefOri_match{1}(x)-prefOri_match{2}(x);
 
end
figure
hist(prefOri_diff);
sgtitle('distribution of change in preferred ori from day 1 to day 2');

%% looking at direction of dfof response in matched cells
resp_avg_match_cell{1}=resp_avg_match{1}(:,:,:,1);
resp_avg_match_cell{1}=squeeze(mean(resp_avg_match_cell{1},1));
resp_avg_match_cell{1}=squeeze(mean(resp_avg_match_cell{1},1));
resp_avg_match_cell{2}=resp_avg_match{2}(:,:,:,1);
resp_avg_match_cell{2}=squeeze(mean(resp_avg_match_cell{2},1));
resp_avg_match_cell{2}=squeeze(mean(resp_avg_match_cell{2},1));
% I now have a 1 by nCell matrix with the average dfof response of each
% cell
figure
subplot(1,2,1)
hist(resp_avg_match_cell{1})
subplot(1,2,2)
hist(resp_avg_match_cell{2})
sgtitle('distribution of average dfof response by cell, day 1 and say 2');
%t-test to determine whether there is a significant difference in average
%response
[h,p]=ttest(resp_avg_match_cell{1},resp_avg_match_cell{2})


% find the timecourses of the supressed cells

% find the contrast response functions of the supressed cells
supressed_subset{1}=con_resp_curves_match{1}(:,find(resp_avg_match_cell{1} < 0))
supressed_subset{2}=con_resp_curves_match{2}(:,find(resp_avg_match_cell{2} < 0))

figure;
subplot(1,2,1)
plot(cons,supressed_subset{1},'LineWidth',.005,'color',[.2 .2 .2]);
xlim([.1 1])
title('Contrast responses 30-degree');
subplot(1,2,2)
plot(cons,supressed_subset{2},'LineWidth',.005,'color',[.2 .2 .2]);
xlim([.1 1])
title('Contrast responses fullfield');
sgtitle('Supressed cells, averaged over orientations');

%% looking at responses for cells with downward contrast resp
down_subset=cell(1,2);
down_subset{1} = resp_avg_match_cell{1}(:,find(con_resp_curves_match{1}(1,:) > con_resp_curves_match{1}(5,:)));
down_subset{2} = resp_avg_match_cell{2}(:,find(con_resp_curves_match{2}(1,:) > con_resp_curves_match{2}(5,:)));


figure;
subplot(1,2,1)
hist(down_subset{1});
title('Contrast responses 30-degree');
subplot(1,2,2)
hist(down_subset{2});
title('Average responses fullfield');
sgtitle('Downward conResp cells, averaged over orientations');


%% comparing matched vs. unmatched cells
%make a set of indices for unmatched cells
inds =  1:size(data_dfof_all{1},3);
x=ismember(inds, match_ind);
inds = inds(~x);

resp_avg_unmatch = resp_avg_all{1}(:,:,inds);
temp_resp_avg_match=resp_avg_match{1}(:,:,:,1);%making a set of just the mean responses for day 1

mean(mean(mean(resp_avg_unmatch)))
mean(mean(mean(temp_resp_avg_match)))
%% naka rushton fit
%fir get data averaged over orientations
con_resp_mat_match=cell(1,nd);
for day = 1:nd
    con_resp_mat_match{day} = squeeze(mean(resp_avg_match{day}(:,:,:,1),1))';
    con_resp_mat_match{day}=max(con_resp_mat_match{day},0);
end

nakaRush_fits = cell(1,nd);
for day = 1:nd
    nakaRush_fits{day} = nakaRushtonFit(con_resp_mat_match{day},cons);
end 
%%
fclose(fid)
%% making an output dataframe from all the response data
outputArray = extract_all_trial_resp(nd,nOri,nCon,nCells,nTrials,resp_cell_match);
%restrict to cells that are responsive on at least one day
for i =1:size(outputArray,1)
    cellID = outputArray(i,1);
    if cellID ~= resp_ind_either;
        outputArray(i,:) = NaN;
    end
end
%removing any extra rows that are all NaN
outputArray = outputArray(all(~isnan(outputArray),2),:); 

%enter cell type, preferred ori and sharpness
for i =1:size(outputArray,1)
    cellID = outputArray(i,1);
    day = outputArray(i,5);
    outputArray(i,6)=red_ind_match(cellID);
    outputArray(i,7)=prefOri_match{day}(cellID);
    outputArray(i,8)=pref_oris_vonM{day}(cellID);
    outputArray(i,9)=vonMises_output{day}(cellID,2);
    outputArray(i,10)=nakaRush_fits{day}(cellID).C50r;
end

%make it into a table
outputTable = array2table(outputArray);
outputTable.Properties.VariableNames(1:10) = {'cellID','dfof','con','ori','day','cell_type','pref_ori','fit_pref_ori','sharpness','C50'};
outputTable.cell_type = categorical( outputTable.cell_type,0:1,{'PYR' 'IN'});
writetable(outputTable,fullfile(fn_multi,'full_dfof_table.csv'));