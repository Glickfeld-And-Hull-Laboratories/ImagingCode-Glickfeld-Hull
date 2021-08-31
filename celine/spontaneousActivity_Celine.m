%will tack this into multiDay_conResp


%% sort TCs into day 1 and day 2, only good cells
%sort into day 1 and day 2, only want red cells
red_TCs = {cellTCs_all{1}(:,red_d1), cellTCs_all{2}(:,red_d2)};


%% spontaneous activity for red cells
% I might adjust this so that it depends on how many red cells there are

pct_events_red = cell(size(red_TCs)); %I don't want pct_events to be overwritten, 
%may need to change this

for iday = 1:size(red_TCs,2) %will loop through each day
    start=1;
    figure
    sgtitle(sprintf(['Red only day ' num2str(iday)]));
    this_day_id = day_id(iday);
    d_all = red_TCs{iday}; %find the data for that day
 
    stimstart = (nOn+1):(nOn+nOff):size(d_all,1)';
    stimon = cell2mat(arrayfun(@(x) x:(x+nOn),stimstart,'unif',0));
    stimoff = setdiff(1:size(d_all,1),stimon);
    d_off = d_all(stimoff,:);
    dff = (d_off-mean(d_off,1))./mean(d_off,1);
    tt = (1:size(d_off,1))./15/60;
    for icell = 1:size(d_all,2)
        if start>5
            cell_lab = icell - 5;
            %print(fullfile(fn_multi,['red_spontaneous_dff_day'  num2str(iday) '_cell_' num2str(cell_lab) 'to' num2str(icell)]),'-dpdf','-fillpage');
            figure;
            sgtitle(sprintf(['Green only day ' num2str(iday)]));
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
    pct_events_red{iday} = sum(dff_test,1)./size(dff,1);
end

%%
nc = cellfun(@(x) size(x,2),pct_events_red);
figure; 
for iday = 1:size(red_TCs,2);
    subplot(1,2,iday); hold on
    d = pct_events_red{iday};
    histogram(d,10);
    vline(mean(d));
    xlabel('pct events');
    ylabel('# cells');
    title(['Spontaneous events day ',num2str(iday)])
end
print(fullfile(fn_multi,'RedSpontaneousEvents'),'-dpdf')

%% spontaneous activity for all cells

pct_events_all = cell(size(cellTCs_all)); 
nexamplecells = 5;

for iday = 1:size(cellTCs_all,2) %will loop through each day
    this_day_id = day_id(iday);
    d_all = cellTCs_all{iday}; %find the data for that day
    
    if size(d_off,2) < nexamplecells
        n = size(d_off,2);
    else
        n = nexamplecells;
    end
    ind = randsample(size(d_all,2),n);

    figure
    sgtitle(sprintf(['Rand subset all cells day ' num2str(iday)]));
    stimstart = (nOn+1):(nOn+nOff):size(d_all,1)';
    stimon = cell2mat(arrayfun(@(x) x:(x+nOn),stimstart,'unif',0));
    stimoff = setdiff(1:size(d_all,1),stimon);
    d_off = d_all(stimoff,:);
    dff = (d_off-mean(d_off,1))./mean(d_off,1);
    tt = (1:size(d_off,1))./expt(this_day_id).frame_rate./60;
    for icell = 1:n
        tc = dff(:,ind(icell));
        subplot(nexamplecells,1,icell)
        plot(tt,tc,'k-','LineWidth',1);
        figXAxis([],'time in expt (min)',[tt(1) tt(end)],0:5:tt(end),0:5:tt(end))
        figYAxis([],'dF/F',[-1 2])
        figAxForm([],0)
        title(sprintf('Cell #%s',num2str(ind(icell))))
    end
   
    dff_3sd = (std(dff) + mean(dff))*3;
    dff_test = dff > dff_3sd;
    pct_events_all{iday} = sum(dff_test,1)./size(dff,1);
    print(fullfile(fn_multi,['all_spontaneous_dff_day'  num2str(iday) '_random_subset']),'-dpdf','-fillpage');
end
%%
nc = cellfun(@(x) size(x,2),pct_events_all);
figure; 
for iday = 1:size(cellTCs_all,2);
    subplot(1,2,iday); hold on
    d = pct_events_all{iday};
    histogram(d,10);
    vline(mean(d));
    xlabel('pct events');
    ylabel('# cells');
    title(['Spontaneous events day ',num2str(iday)])
end
print(fullfile(fn_multi,'AllSpontaneousEvents'),'-dpdf')