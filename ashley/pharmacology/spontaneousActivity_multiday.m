%% which dataset and parameters
clear all
clear global
close all
ds = 'DART_V1_PV_contrast'; %dataset info
dataStructLabels = 'contrastxori';

rc = behavConstsAV; %directories
eval(ds)


%% which days to use?
day1_id = 8;
day2_id = 9;
day3_id = nan;

%% mouse and saving info
mouse = expt(day1_id).mouse;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging');
fnout = fullfile(fn,'multiday');
mkdir(fnout)
%% load multiday cell TCs
load(fullfile(fnout,'timecourses.mat'))
%% sort TCs into day 1 and day 2, only good cells
alignedCells_tag = cell2mat({red_cellTCs.pass})==1;
alignedCells_notag = cell2mat({green_cellTCs.pass})==1;

data_tag = cat(2,{cell2mat({red_cellTCs.tc_day1})},{cell2mat({red_cellTCs.tc_day2})});
data_notag = cat(2,{cell2mat({green_cellTCs.tc_day1})},{cell2mat({green_cellTCs.tc_day2})});

% day 1
expDate_day1 = expt(day1_id).date;
expTime_day1 = eval(['expt(day1_id).' dataStructLabels '_time']);

% day 2
expDate_day2 = expt(day2_id).date;
expTime_day2 = eval(['expt(day2_id).' dataStructLabels '_time']);

if ~isnan(day3_id)
    % day 3
    expDate_day3 = expt(day3_id).date;
    expTime_day3 = eval(['expt(day3_id).' dataStructLabels '_time']);
end
%% load mworks data
mw_day1 = loadMworksFile(mouse,expDate_day1,expTime_day1{1});
mw_day2 = loadMworksFile(mouse,expDate_day2,expTime_day2{1});
if ~isnan(day3_id)
    mw_day3 = loadMworksFile(mouse,expDate_day3,expTime_day3{1});
    mw = {mw_day1,mw_day2,mw_day3};
else
    mw = {mw_day1,mw_day2};
end

%% 
% data = data_tag;

%% spontaneous activity
nexamplecells = 5;
pct_events = cell(2,size(data_tag,2));
tag_id = {'tag','notag'};
for itag = 1:2
    d_all = eval(['data_' tag_id{itag}]);
    if size(d_off,2) < nexamplecells
        n = size(d_off,2);
    else
        n = nexamplecells;
    end
    ind = randsample(size(d_all{1},2),n);

    for iday = 1:size(d_all,2)
        dt = eval(['expt(day' num2str(iday),'_id).date']);
        figure
        suptitle([dt ', ' tag_id{itag} ', random cells'])
        d = d_all{iday};
        mw = eval(['mw_day' num2str(iday)]);
        on = mw.nScansOn;
        off = mw.nScansOff;
        stimstart = (on+1):(on+off):size(d,1)';
        stimon = cell2mat(arrayfun(@(x) x:(x+on),stimstart,'unif',0));
        stimoff = setdiff(1:size(d,1),stimon);
        d_off = d(stimoff,:);
        dff = (d_off-mean(d_off,1))./mean(d_off,1);
        tt = (1:size(d_off,1))./expt(day1_id).frame_rate./60;
        for icell = 1:n
            tc = dff(:,ind(icell));
            subplot(nexamplecells,1,icell)
            plot(tt,tc,'k-','LineWidth',1);
            figXAxis([],'time in expt (min)',[tt(1) tt(end)],0:5:tt(end),0:5:tt(end))
            figYAxis([],'dF/F',[-1 2])
            figAxForm([],0)
            title(sprintf('Cell #%s',num2str(ind(icell))))
        end
        print(fullfile(fnout,['spontaneous_dff_' dt '_' tag_id{itag}]),'-dpdf','-fillpage')
        dff_3sd = (std(dff) + mean(dff))*3;
        dff_test = dff > dff_3sd;
        pct_events{itag,iday} = sum(dff_test,1)./size(dff,1);
    end
end
nc = cellfun(@(x) size(x,2),pct_events);
figure; 
for itag = 1:2
    subplot(1,2,itag); hold on
    for iday = 1:size(data_tag,2)
        d = pct_events{itag,iday};
        plot(ones(nc(itag,iday),1)*iday,d,'k.','MarkerSize',10)
        errorbar(iday,mean(d),ste(d,2),'b.','MarkerSize',20)
    end
    x_labels = cell(size(data_tag));
    for i = 1:size(data_tag,2)
        ind = expt(day1_id).multiday_matchdays(i);
        x_labels{i} = [num2str(expt(ind).multiday_timesincedrug_hours) 'hrs'];
    end
    figXAxis([],'',[0 4],1:size(data_tag,2),x_labels)
    figYAxis([],'spontaneous events',[])
    figAxForm
    title(tag_id{itag})
end
print(fullfile(fnout,'spontaneousEvents'),'-dpdf')