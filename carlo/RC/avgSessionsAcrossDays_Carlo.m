%expt info
mouse = '1081';
d = cell(1,2);
d{1}= '210214';
d{2} = '210215';
cue = ['0','1'];
cmat = ['r', 'k'];
base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\carlo\analysis\2P';
frameRateHz = 30;
prewin_frames = round(1500./frameRateHz);
postwin_frames = round(3000./frameRateHz);
tt = (-prewin_frames:postwin_frames-1).*1000/frameRateHz;

%load data for each cue condition on each day
figure;
for ic = 1:length(cue)
    c = cue(ic);
    data = [];
    for id = 1:length(d)
        date = d{id};
        expt = [date c];
        if id==2 & ic==1 || id==1 & ic==2
            t = 'getTC_001';
        else
            t = 'getTC_000';
        end
        fn = fullfile(base,[expt '_img' mouse], t,['figDataUntrimmed_' expt '_img' mouse '.mat']);
        data_temp = load(fn);
        %average across trials and convert to Hz
        data_avg = mean(data_temp.targetAlign_events,3)./(1./frameRateHz);
        %concatenate cells from each day
        data = [data data_avg];
        size(data)
    end
    shadedErrorBar(tt, mean(data,2), std(data,[],2)./sqrt(size(data,2)),cmat(ic));
    hold on
end
vline(0)
vline(770)
xlim([-1000 2000])
title([mouse ' dates: ' d{1} ' ' d{2}])
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\public\for Court\' mouse '_blockedCueSession.pdf'],'-dpdf')
savefig(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\public\for Court\' mouse '_blockedCueSession.fig'])

        