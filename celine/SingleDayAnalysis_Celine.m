 % to perform additional analyses on already extracted timecourses
clear all; clear global; close all; clc
%ds = 'con_ori_nonDART';
ds = 'DART_V1_contrast_ori_Celine'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds)


day_id = 280;
% identifying animal and run
mouse = expt(day_id).mouse;
date = expt(day_id).date;

imgFolder = expt(day_id).contrastxori_runs{1};

frame_rate = 15;

%setting my paths
fn_base = 'Z:\home\ACh\Analysis\2p_Analysis\';
fn = fullfile(fn_base,mouse,date,imgFolder);
cd(fn);

load('input.mat')
load('mask_cell.mat')
load('regOuts&Img.mat')
load('TCs.mat')
% get info
nCells = size(npSub_tc,2);
nOn = input.nScansOn;
nOff = input.nScansOff;
ntrials = size(input.tGratingContrast,2);
tCon = celleqel2mat_padded(input.tGratingContrast);
cons = unique(tCon);
nCon = length(cons);
ind_con = find(tCon == max(cons(:)));
tDir = celleqel2mat_padded(input.tGratingDirectionDeg);
tOri = tDir;
tOri(find(tDir>=180)) = tDir(find(tDir>=180))-180;
oris = unique(tOri);
nOri = length(oris);

%%

% % % plot cells with ID numbers to see if there are any you want to get ride of
% cell_stats=regionprops(mask_cell);
% figure; imagesc(mask_cell)
% hold on
% bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
% plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',.5); hold on;
% for iC = 1:length(find(mask_cell))
%     text(cell_stats(iC).Centroid(1), cell_stats(iC).Centroid(2), num2str(iC), 'Color', 'white',...
%             'Fontsize', 10, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
%     
% end
% % remove cells you don't want
% ind=5;
% mask_cell(find(mask_cell == ind))=0;
% mask_cell_red(find(mask_cell_red == ind))=0;
% 
% data_tc(:,ind)=[];
% mask_label(:,ind)=[];
% npSub_tc(:,ind)=[];
% nCells=nCells-1;
%get average FOV with masks

% 
figure;
imagesc(data_avg);
colormap gray
title('average FOV');
hold on
bound = cell2mat(bwboundaries(mask_cell(:,:,1)));
plot(bound(:,2),bound(:,1),'.','color','g','MarkerSize',2);
bound = cell2mat(bwboundaries(mask_cell_red(:,:,1)));
%plot(bound(:,2),bound(:,1),'.','color','r','MarkerSize',2);
hold off
print(fullfile(fn,['FOV_with_masks']),'-dpdf');

% reshape by trials and look at responses
%getting df/f for each trial, using a baseline window
%% new way
stimStart = (nOff/2)+1; %this indicates both the perdiod to trim off the start and the stim on period after trimming
stimEnd=stimStart+nOn-1;

resp_win = stimStart:stimEnd;
base_win = 1: (nOff/2);

npSub_tc_new=npSub_tc(stimStart:size(npSub_tc,1),:); %trimming npSub_tc to correct size( sacrifiace one trial)
npSub_tc_new=padarray(npSub_tc_new,30,0,'post');
data_tc_trial = reshape(npSub_tc_new, [nOn+nOff,ntrials,nCells]);
data_f_trial = mean(data_tc_trial(1:30,:,:),1);
data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);


%% split into trials
data_resp = zeros(nCells, nOri, nCon,2);
h = zeros(nCells, nOri, nCon);
p = zeros(nCells, nOri, nCon);
tCon = tCon(:,1:ntrials);
tOri = tOri(:,1:ntrials);
rect_dfof_trial=data_dfof_trial;
rect_dfof_trial(find(rect_dfof_trial<0))=0;

for iOri = 1:nOri
    ind_ori = find(tOri == oris(iOri));
    for iCon = 1:nCon
        ind_con = find(tCon == cons(iCon));
        ind = intersect(ind_ori,ind_con); %for every orientation and then every contrast, find trials with that con/ori combination
        data_resp(:,iOri,iCon,1) = squeeze(mean(mean(rect_dfof_trial(resp_win,ind,:),1),2));
        data_resp(:,iOri,iCon,2) = squeeze(std(mean(rect_dfof_trial(resp_win,ind,:),1),[],2)./sqrt(length(ind)));
        [h(:,iOri,iCon), p(:,iOri,iCon)] = ttest(mean(rect_dfof_trial(resp_win,ind,:),1), mean(rect_dfof_trial(base_win,ind,:),1),'dim',2,'tail','right','alpha',0.05./(nOri.*3-1));
    end
end

h_all = sum(sum(h,2),3);

resp=logical(h_all);
%
%pref_ori = zeros(1,nCells);
%orth_ori = zeros(1,nCells);
pref_con = zeros(1,nCells);
data_ori_resp = zeros(nCells,nOri); %at pref con
data_con_resp = zeros(nCells,nCon); %at pref ori
data_orth_resp=zeros(nCells,1);

%I want to pull out the responses for each cell at it's preferred orientations, for
%all contrasts, and at it's preferred contrast, for all orientations
 [max_val, pref_ori] = max(mean(data_resp(:,:,:,1),3),[],2);
 for iCell = 1:nCells
      
      [max_val_con, pref_con(1,iCell)] = max(squeeze(mean(data_resp(iCell,:,:,1),2))',[],2);
      if pref_ori(iCell)<= 4
          orth_ori(iCell)=pref_ori(iCell)+4;
      elseif pref_ori(iCell)>4
          orth_ori(iCell)=pref_ori(iCell)-4;
      end
      data_ori_resp(iCell,:)=data_resp(iCell,:,pref_con(iCell),1);
%      data_orth_resp(iCell,:)=mean(data_resp(iCell,orth_ori(iCell),:,1));
      data_con_resp(iCell,:)=data_resp(iCell,pref_ori(iCell),:,1);
end

% %% calculating OSI 
% %OSI only seems to work on rectified dfof
% OSI1 = (max_val-data_orth_resp)./(max_val+data_orth_resp);
% movegui('center')
% hist(OSI1)

%% get basic counts
green_inds = 1:nCells;
green_inds = setdiff(green_inds, find(mask_label));


nGreen = length(green_inds);
GreenResp = intersect(green_inds,find(resp));
nGreenResp = length(GreenResp);

nRed = length(find(mask_label));
RedAll =find(mask_label);%to take all red cells not only the responsive ones
RedResp =intersect(RedAll,find(resp));
nRedResp = length(RedResp);

% comparing the number of stimuli that each cell was responsive to
%for only the responsive cells, how many sitmulus condiditions did they
%respond to?
PC_nResp = h_all(GreenResp);
red_nResp = h_all(RedAll);

%isolate dfof values for conditions each cell was responsive to
dfof_resp = data_resp(:,:,:,1) .* h;
dfof_resp(dfof_resp==0) = NaN;
dfof_resp_green = nanmean(nanmean(dfof_resp(GreenResp,:,:),3),2);
dfof_resp_red = nanmean(nanmean(data_resp(RedAll,:,:,1),3),2);


%how many different orientations did each cell respond to, regardless of
%contrast?
PC_nOris_all = sum(logical(sum(h(GreenResp,:,:),3)),2);
red_nOris_all = sum(logical(sum(h(RedAll,:,:),3)),2);


PC_nCons_all = sum(logical(squeeze(sum(h(GreenResp,:,:),2))),2);
red_nCons_all = sum(logical(squeeze(sum(h(RedAll,:,:),2))),2);



% make table

countsTable = table([nGreen; nGreenResp;mean(PC_nResp);mean(PC_nOris_all);mean(PC_nCons_all);mean(dfof_resp_green)],[nRed;nRedResp;mean(red_nResp);mean(red_nOris_all);mean(red_nCons_all);mean(dfof_resp_red)],'VariableNames',{'PCs' 'INs'}, 'RowNames',{'Total cells' 'Responsive cells','Avrg conditions resp','Avrg oris resp','Avrg cons resp','Avrg dfof of resp'})
writetable(countsTable,fullfile(fn,'counts.csv'),'WriteRowNames',true)
% looking at time courses for all cells (responsive and non-responsive)
% data_f_trial = mean(data_tc_trial(nOff/2:nOff,:,:),1);
% data_dfof_trial = bsxfun(@rdivide, bsxfun(@minus,data_tc_trial, data_f_trial), data_f_trial);

%% print histogram of preferred contrasts
% pref_con_converted = cons(pref_con);
% figure
% subplot(2,1,1)
% hist(pref_con_converted(:,GreenResp))
% subplot(2,1,2)
% hist(pref_con_converted(:,RedAll))
% print(fullfile(fn,['pref_con_hist']),'-dpdf');
%% narrow down to the preferred ori and preferred contrast for each cell
tc_trial_avrg=cell(1,nCon);

for iCon = 1:nCon

    tc_trial_avrg_temp=nan(nOn+nOff,nCells);
    for i=1:nCells
        temp_TCs=data_dfof_trial(:,:,i);
        %identify the trials where ori = pref ori
        temp_ori = oris(pref_ori(i)); %find the preferred ori of this cell and convert to degrees
        ori_inds = find(tDir==temp_ori); %these are the trials at that ori
        %identify which contrasts the cell responded to
        %temp_h = squeeze(h(i,pref_ori(i),:)); %pull out the h values for the preferred ori and this cell only
        %temp_cons = (temp_h').*cons;
        %temp_cons(temp_cons==0)=[];
        con_inds=find(tCon==cons(iCon));
        temp_trials = intersect(ori_inds, con_inds);
        tc_trial_avrg_temp(:,i)=mean(temp_TCs(:,temp_trials),2);
        
    end
tc_trial_avrg{iCon}=tc_trial_avrg_temp;
end

%% plotting
%tc_trial_avrg = squeeze(mean(data_dfof_trial,2));%average over trials, one row per cell

for iCon = 1:nCon
    tc_green = tc_trial_avrg{iCon}(:,GreenResp);
    tc_red = tc_trial_avrg{iCon}(:,RedResp);
    tc_green_avrg=cell(1,2);
    tc_red_avrg=cell(1,2);
    tc_green_avrg{1} = mean(tc_green,2);%average tc for responsive green cells
    tc_green_avrg{2}=std(tc_green,[],2);
    tc_red_avrg{1} = mean(tc_red,2); %average tc for all red cells
    tc_red_avrg{2}=std(tc_red,[],2);
    
    %convert to se 
    tc_green_avrg{2} = tc_green_avrg{2}/sqrt(size(GreenResp,1));
    tc_red_avrg{2} = tc_red_avrg{2}/sqrt(size(RedResp,1));
    
    
    
    figure
    x=1:(size(tc_green,1));
    x=(x-30)/15;
    shadedErrorBar(x,tc_red_avrg{1},tc_red_avrg{2},'r');
    %plot(x,tc_red,'r');
    ylim([-.05 .3]);
    hold on
    shadedErrorBar(x,tc_green_avrg{1},tc_green_avrg{2});
    title(num2str(cons(iCon)));
    saveas(gcf,fullfile(fn,['contrast_',num2str(cons(iCon)), '_tcPlot.jpg']));

end


%% plotting the contrasts together
figure;
for iCon = 1:nCon
    tc_green = tc_trial_avrg{iCon}(:,GreenResp);
    tc_red = tc_trial_avrg{iCon}(:,RedResp);
    tc_green_avrg=cell(1,2);
    tc_red_avrg=cell(1,2);
    tc_green_avrg{1} = mean(tc_green,2);%average tc for responsive green cells
    tc_green_avrg{2}=std(tc_green,[],2);
    tc_red_avrg{1} = mean(tc_red,2); %average tc for all red cells
    tc_red_avrg{2}=std(tc_red,[],2);
    
    %convert to se 
    tc_green_avrg{2} = tc_green_avrg{2}/sqrt(size(GreenResp,1));
    tc_red_avrg{2} = tc_red_avrg{2}/sqrt(size(RedResp,1));
    

    
   
    x=1:(size(tc_green,1));
    x=(x-30)/15;
    if iCon == 1
        shadedErrorBar(x,tc_red_avrg{1},tc_red_avrg{2},'b');
        %shadedErrorBar(x,tc_green_avrg{1},tc_green_avrg{2},'b');
        hold on
    elseif iCon== 2
        shadedErrorBar(x,tc_red_avrg{1},tc_red_avrg{2},'g');
        %shadedErrorBar(x,tc_green_avrg{1},tc_green_avrg{2},'g');
        hold on
    elseif iCon==3
        shadedErrorBar(x,tc_red_avrg{1},tc_red_avrg{2},'c');
        %shadedErrorBar(x,tc_red_avrg{1},tc_green_avrg{2},'c');
        hold on
    end
    %plot(x,tc_red,'r');
    ylim([-.05 .3]);
    
    
    %shadedErrorBar(x,tc_green_avrg{1},tc_green_avrg{2});
    
    %saveas(gcf,fullfile(fn,['contrast_',num2str(cons(iCon)), '_tcPlot.jpg']));

end

%% to see what the red cells look like individually
figure;
plot(tc_red);
title('Timecourses for all red cells stimuli');
legend
%% raw timecourses for all stimuli / all cells

tc_cell_avrg_all = mean(data_dfof_trial,3);%average pver cells, one row per trial
tc_trial_avrg_all = squeeze(mean(data_dfof_trial,2));%average over trials, one row per cell
tc_cell_trial_avrg_all = mean(tc_cell_avrg_all,2);%average over trials and cells

figure;
plot(tc_trial_avrg_all(:,green_inds), 'LineWidth',.005,'color',[.25 .25 .25]);
hold on;
plot(tc_trial_avrg_all(:,RedAll), 'LineWidth',.005, 'color','r');
hold on
plot(tc_cell_trial_avrg_all, 'LineWidth',2, 'color','k');
title('Timecourses all cells all stimuli');
vline(30,'g')
ylim([-.1 .3]);
hold off
saveas(gcf,fullfile(fn,[mouse '-' date 'raw_tc.jpg']));

%% plotting contrast and ori functions 
keepCells=union(GreenResp,RedAll);
nKeep=size(keepCells,1);
[RedAll,red_ind_keep] = intersect(keepCells,RedAll,'stable');

if nKeep<36
    [n n2] = subplotn(nKeep);
    tot = n.*n2;
else
    n = 6;
    n2 = 6;
    tot = 36;
end

% plot ori tuning at each contrast
figure;
movegui('center')
start = 1;

for iCell = 1:nKeep
    if start>tot
        figure; movegui('center')
        start = 1;
    end
    subplot(n,n2,start)


    if find(find(iCell==red_ind_keep))
          for iCon = 1:nCon
            errorbar(oris, data_resp(iCell,:,iCon,1), data_resp(iCell,:,iCon,2),'-r')
            hold on
            title (['Cell ' num2str(iCell)]);
          end
        
    else
         for iCon = 1:nCon
            errorbar(oris, data_resp(iCell,:,iCon,1), data_resp(iCell,:,iCon,2),'-k')
            hold on
            title (['Cell ' num2str(iCell)]);
         end
    end
    start= start+1;
    %ylim([-0.1 inf])

   
end


% % 
% % % plots contrast preference at preferred orientation
% % figure;
% % movegui('center')
% % start = 1;
% % for iCell = 1:nCells
% %     if start>tot
% %         figure; movegui('center')
% %         start = 1;
% %     end
% %     subplot(n,n2,start)
% % 
% %     errorbar(cons, squeeze(data_resp(iCell,pref_ori(iCell),:,1)), squeeze(data_resp(iCell,pref_ori(iCell),:,2)),'-o')
% %     ylim([-0.1 inf])
% %     start = start+1;
% % 
% % end
% % 
% % %% naka Rushton fit at pref ori
% % nakaRush_fits = nakaRushtonFit(data_con_resp,cons);
% % x=squeeze(struct2cell(nakaRush_fits));
% % 
% % for i = 1:nCells
% %    c50(i)=x{3,i}; 
% % end
% % c50_PCresp=c50(GreenResp)';
% % c50_INresp=c50(RedAll)';
% %save(fullfile(fn,'nakaRush_fits.mat'),'nakaRush_fits');
% %save(fullfile(fn,'nakaRush_fits_PCresp.mat'),'nakaRush_fits_PCresp');
% %save(fullfile(fn,'nakaRush_fits_ICresp.mat'),'nakaRush_fits_ICresp');
% 
% %% go back to my own von M fit?
% 
% % %1 initial von M fit on real data/real orientations collected
% % vonMises_output = nan(nCells, 6);
% % 
% % for i = 1:nCells
% %     thisCell = data_ori_resp(i,:);
% %     [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(oris,thisCell);
% %     vonMises_output(i,:)=[b_hat, k1_hat, R1_hat,u1_hat,sse,R_square];
% % end
% % 
% % clear thisCell b_hat k1_hat R1_hat u1_hat sse R_square %get rid of the last iteration of these looping variables
% % save(fullfile(fn,'vonM_output.mat'),'vonMises_output')
% % 
% % PC_k = vonMises_output(GreenResp,2);
% % IN_k = vonMises_output(RedAll,2);
% 
% % save(fullfile(fn,'vonM_PCresp.mat'),'vonM_PCresp')
% % save(fullfile(fn,'vonM_INresp.mat'),'vonM_INresp')
% 
% 
% %2 fit with higher resolution
% 
% %3 bootstrap to find reliability
% %for each cell
%     %for each boot
%         %redo finding the average response to each orientation,still only at the preferred contrast. Feed that new dataset into 
%         %find the fit from raw data
%         %use those values to find the fit at higher resolution
%         %find the preferred ori from the step above
%         %comapre the preferred ori to the original preferred ori from 2
%         %above
%      %end
%      %sort the difference value low to high
%      %find the one in the 900th position - save this
% %end
% 
% %determine which cells have a value less than 22.5, ie are reliable
% % %% making overall output table
% % PC_responses = [PC_nResp PC_nOris_all PC_nCons_all dfof_resp_green c50_PCresp PC_k];
% % IN_responses = [red_nResp red_nOris_all red_nCons_all dfof_resp_red c50_INresp IN_k];
% % 
% % PC_table=array2table(PC_responses, 'VariableNames',{'nResp' 'nOris' 'nCons' 'dfof_resp' 'c50' 'k'});
% % IN_table=array2table(IN_responses, 'VariableNames',{'nResp' 'nOris' 'nCons' 'dfof_resp' 'c50' 'k'});
% % 
% % writetable(PC_table,fullfile(fn,'PC_responses.csv'));
% % writetable(IN_table,fullfile(fn,'IN_responses.csv'));
% 
% 
% %%
% 
% % %% trying LG fit
% % %for this I need the timecourses at the peferred contrast
% % [avgResponseEaOri,semResponseEaOri,vonMisesFitAllCellsAllBoots,fitReliability,R_square,tuningTC,sharpness]...
% %     = getOriTuningLG_CC(npSub_tc,input,5);
% % %% extracting output from vonM fits
% % % the output is in a 181 orientations X 1001 bootstraps X nCells matrix,
% % % and there is a fitReliability output that tells the 90th percentile of
% % % the difference between the bootstrapped fits and the original fit.
% % % Ultimately I probably only want to take cells that have a fitReliability
% % % value below 22.5 but right now none of them do
% % median(fitReliability,2)
% % 
% % prefOris = nan(size(vonMisesFitAllCellsAllBoots,2),nCells);
% % ogFitPref = nan(1,nCells);
% % for i = 1:nCells
% %     for fit = 2:size(vonMisesFitAllCellsAllBoots,2)
% %         prefOris(fit,i) = mean(find(vonMisesFitAllCellsAllBoots(:,fit,i)==max(vonMisesFitAllCellsAllBoots(:,fit,i))));
% %         ogFitPref(i) = mean(find(vonMisesFitAllCellsAllBoots(:,1,i)==max(vonMisesFitAllCellsAllBoots(:,1,i))));
% %     end
% % end
% % meanPref = mean(prefOris,1);
% % %I now have meanPref which tells me the average orientation indicated by
% % %the 1000 bootstraps, and the max from the original fit. I'm not sure how
% % %to get the oringal U value, or the K value
% % vonMisesFitAllCells = squeeze(vonMisesFitAllCellsAllBoots(:,1,:)); % takes just the first of the bootstraps, which is the original fit 
% % [maxResp prefOri] = max(vonMisesFitAllCells,[],1);
% % 
% % %combine max Resp, prefOri, and fit reliability 
% % vonM_output = [maxResp; prefOri;fitReliability];
% % 
% % save(fullfile(fn,'vonM_fits.mat'),'vonMisesFitAllCells');
% % save(fullfile(fn,'vonM_output.mat'),'vonM_output');
% % 
% 

%% look at spontaneous actiivty - HT+ cells
RedAll =find(mask_label);
red_TCs = npSub_tc(:,RedAll);
start=1;
figure
sgtitle(sprintf(['Red only']));
d_all = red_TCs; %find the data for that day

pct_events_red=cell(1,1);
stimstart = (nOn+1):(nOn+nOff):size(d_all,1)';
stimon = cell2mat(arrayfun(@(x) x:(x+nOn),stimstart,'unif',0));
stimoff = setdiff(1:size(d_all,1),stimon);
d_off = d_all(stimoff,:);
dff = (d_off-mean(d_off,1))./mean(d_off,1);
tt = (1:size(d_off,1))./frame_rate./60;
for icell = 1:size(d_all,2)
    if start>5
        cell_lab = icell - 5;
        print(fullfile(fn,['red_spontaneous_dff_cell_' num2str(cell_lab) 'to' num2str(icell)]),'-dpdf','-fillpage');
        figure;
        sgtitle(sprintf(['Red only']));
        movegui('center')
        start = 1;
    end
    tc = dff(:,icell);
    subplot(5,1,start)
    plot(tt,tc,'k-','LineWidth',1);
    figXAxis([],'time in expt (min)',[tt(1) tt(end)],0:5:tt(end),0:5:tt(end))
    figYAxis([],'dF/F',[-1 2])
    figAxForm([],0)
    title(sprintf('Cell #%s',num2str(icell)))
    start = start+1;
end

dff_3sd = (std(dff) + mean(dff))*3;
dff_test = dff > dff_3sd;
pct_events_red{1}= sum(dff_test,1)./size(dff,1);


nc = cellfun(@(x) size(x,2),pct_events_red);
figure; 
d = pct_events_red{1};
histogram(d,10);
vline(mean(d));
xlabel('pct events');
ylabel('# cells');
title(['Spontaneous events'])
print(fullfile(fn,'RedSpontaneousEvents'),'-dpdf')

% spontaneous events - HT-
green_inds = 1:nCells;
green_inds = setdiff(green_inds, find(mask_label));
green_TCs = npSub_tc(:,green_inds);

start=1;
figure
sgtitle(sprintf(['Green only']));
d_all = green_TCs; %find the data for that day

pct_events_green=cell(1,1);
stimstart = (nOn+1):(nOn+nOff):size(d_all,1)';
stimon = cell2mat(arrayfun(@(x) x:(x+nOn),stimstart,'unif',0));
stimoff = setdiff(1:size(d_all,1),stimon);
d_off = d_all(stimoff,:);
dff = (d_off-mean(d_off,1))./mean(d_off,1);
tt = (1:size(d_off,1))./frame_rate./60;
% for icell = 1:size(d_all,2)
%     if start>5
%         cell_lab = icell - 5;
%         print(fullfile(fn,['green_spontaneous_dff_cell_' num2str(cell_lab) 'to' num2str(icell)]),'-dpdf','-fillpage');
%         figure;
%         sgtitle(sprintf(['Green only']));
%         movegui('center')
%         start = 1;
%     end
%     tc = dff(:,icell);
%     subplot(5,1,start)
%     plot(tt,tc,'k-','LineWidth',1);
%     figXAxis([],'time in expt (min)',[tt(1) tt(end)],0:5:tt(end),0:5:tt(end))
%     figYAxis([],'dF/F',[-1 2])
%     figAxForm([],0)
%     title(sprintf('Cell #%s',num2str(icell)))
%     start = start+1;
% end

dff_3sd = (std(dff) + mean(dff))*3;
dff_test = dff > dff_3sd;
pct_events_green{1}= sum(dff_test,1)./size(dff,1);


nc = cellfun(@(x) size(x,2),pct_events_green);
figure; 
d = pct_events_green{1};
histogram(d,10);
vline(mean(d));
xlabel('pct events');
ylabel('# cells');
title(['Spontaneous events'])
print(fullfile(fn,'GreenSpontaneousEvents'),'-dpdf')




