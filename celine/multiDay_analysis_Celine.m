cd(fn_multi)

PC_ind_match = find(red_ind_match ==0);

%%spontaneous activity
%was there a significant different in spontaneous (not visually evoked)
%activity between days 1 and 2

%isolate the non-HT cells
PC_match_fractActive = cell(1,2);
for i=1:2
   PC_match_fractActive{i} = fractTimeActive_match{i}(PC_ind_match);
   mean(PC_match_fractActive{i});
end

[h p] = ttest(PC_match_fractActive{1},PC_match_fractActive{2})
mean(PC_match_fractActive{1})
mean(PC_match_fractActive{2})


%isolate the red matched cells
red_match_fractActive = cell(1,2);
for i=1:2
   red_match_fractActive{i} = fractTimeActive_match{i}(red_ind_match);
   mean(red_match_fractActive{i});
end

[h p] = ttest(red_match_fractActive{1},red_match_fractActive{2})
mean(red_match_fractActive{1})
mean(red_match_fractActive{2})

%% comparing the number of stimuli that each cell was responsive to
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


%% get preferred orientation (and sharpness) from LG fit
%this will produce 1000 bootstraps

avgResponseEaOri_match = cell(1,nd);
semResponseEaOri_match = cell(1,nd);
vonMisesFitAllCellsAllBoots_match = cell(1,nd);
fitReliability_match = cell(1,nd);
R_square_match = cell(1,nd);
tuningTC_match = cell(1,nd);
sharpness_match = cell(1,nd);

for day = 1:nd
    [avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,...
        fitReliability,R_square,tuningTC,sharpness] = getOriTuningLG_CC(cellTCs_match{day}(:,resp_ind_either),input(:,day),5);
    avgResponseEaOri_match{day} = avgResponseEaOri;
    semResponseEaOri_match{day} = semResponseEaOri;
    vonMisesFitAllCellsAllBoots_match{day} = vonMisesFitAllCellsAllBoots;
    fitReliability_match{day} = fitReliability;
    R_square_match{day} = R_square;
    tuningTC_match{day} = tuningTC;
    sharpness_match{day} = sharpness;
end
%% %% extracting output from vonM fits
nRespCell = size(resp_ind_either,1);
%look at median fit reliability for each day
nanmedian(fitReliability_match{1},2)
nanmedian(fitReliability_match{2},2)

for day = 1:nd
    prefOris{day} = nan(size(vonMisesFitAllCellsAllBoots_match{day},2),nRespCell);
    ogFitPref{day} = nan(1,nRespCell);
    for i = 1:nRespCell
        for fit = 2:size(vonMisesFitAllCellsAllBoots_match{day},2)
            prefOris{day}(fit,i) = mean(find(vonMisesFitAllCellsAllBoots_match{day}(:,fit,i)==max(vonMisesFitAllCellsAllBoots_match{day}(:,fit,i))));
            ogFitPref{day}(i) = mean(find(vonMisesFitAllCellsAllBoots_match{day}(:,1,i)==max(vonMisesFitAllCellsAllBoots_match{day}(:,1,i))));
        end
    end
    meanPref{day} = mean(prefOris{day},1);
    %I now have meanPref which tells me the average orientation indicated by
    %the 1000 bootstraps, and the max from the original fit. I'm not sure how
    %to get the oringal U value, or the K value
    vonMisesFitAllCells{day} = squeeze(vonMisesFitAllCellsAllBoots_match{day}(:,1,:)); % takes just the first of the bootstraps, which is the original fit 
    [maxResp{day} prefOri{day}] = max(vonMisesFitAllCells{day},[],1);

    %combine max Resp, prefOri, and fit reliability 
    vonM_output{day} = [maxResp{day}; prefOri{day};fitReliability_match{day}];
end
% 
% save(fullfile(fn,'vonM_fits.mat'),'vonMisesFitAllCells');
% save(fullfile(fn,'vonM_output.mat'),'vonM_output');


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

%% naka rushton fit
%averaged over orientations
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