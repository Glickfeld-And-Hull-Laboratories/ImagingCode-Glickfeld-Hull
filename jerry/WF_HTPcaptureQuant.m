clear all
close all
mice = {'i3312'};
ratios = getCaptureValues_annulus_jerry(mice);

%% plots

nMice = size(ratios,1);

for iMouse = 1:nMice
    currData = ratios{iMouse,2};
    figure();
    plot(currData(:,2),currData(:,1),'-o');
    title(ratios{iMouse,1});
    ylim([0.8 2]);
    Ax = gca;
    Ax.XTick = currData(:,2);
    xlabel('Time since Infusion (hr)');
    ylabel('Capture Index (AU)');
end