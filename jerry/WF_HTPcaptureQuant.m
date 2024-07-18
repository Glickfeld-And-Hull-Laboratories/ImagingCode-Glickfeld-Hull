clear all
close all
mice = {'i3312'};
dat = getCaptureValues_annulus_jerry(mice);
nMice = size(dat,1);
%% plot ratio


for iMouse = 1:nMice
    currData = dat{iMouse,2};
    figure();
    plot(currData(:,2),currData(:,1),'-o');
    title([dat{iMouse,1} ' Ratio']);
    ylim([0.8 2]);
    Ax = gca;
    Ax.XTick = currData(:,2);
    xlabel('Time since Infusion (hr)');
    ylabel('Capture Index (AU)');
end

%% plot unnormalized intensity


for iMouse = 1:nMice
    currData = dat{iMouse,2};
    figure();
    plot(currData(:,2),currData(:,3),'-o');
    title([dat{iMouse,1} ' HTP ROI']);
    ylim([70 200])
    Ax = gca;
    Ax.XTick = currData(:,2);
    xlabel('Time since Infusion (hr)');
    ylabel('Mean HTP ROI Intensity (AU)');
end
