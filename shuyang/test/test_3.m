MWorks_output = load([MWorks_source 'data-i' days '.mat' ]);
MWorks_output = MWorks_output.input;
dfOvF_temp = load([image_dest,'_dfOvF.mat']);
dfOvF = dfOvF_temp.dfOvF_btmbase;
%dfOvF = dfOvF_temp.dfOvF_avglaseroff;
% laseroff = 300;
% laseron = 50;
laseroff = double(MWorks_output.nScansOff);
laseron = double(MWorks_output.nScansOn);
%which frames are the first frame when laser is turned on each time: 151,311...
%ntrials = double(MWorks_output.trialSinceReset);
ntrials = 15;
laseron_frame = zeros(1,ntrials-1);
for i = 1:ntrials-1
    laseron_frame(i) = (laseroff+laseron)*(i-1)+1+laseroff; % the indexes of laser onset, not the whole time of laseron
end

laser_align_dfOvF = zeros(ntrials-1,10+laseron+15,size(dfOvF,1));% trial*frame*ROI, 1s before laser on and 1.5s after laser on 
for r = 1:size(dfOvF,1)
    for f = 1:length(laseron_frame)
        laser_align_dfOvF(f,:,r) = dfOvF(r,laseron_frame(f)-9:laseron_frame(f)+laseron+15);
    end
end

align_fig = figure; 
align_fig.Units = 'centimeters';
align_fig.Position = [1 1 12 3.45];
x = double(1:(10+laseron+15));
x2 = (x./10)-0.9;
colorcode = {'r','d'};
for r = 1:size(dfOvF,1)
    shadedErrorBar(x2,mean(laser_align_dfOvF(:,:,r),1),std(laser_align_dfOvF(:,:,r),0,1)/sqrt(size(laser_align_dfOvF(:,:,r),1)));hold on; 
end
xlabel('time from laser on(s)');
ylabel('df/F');
title(sessions);
ylim([0 0.5]);
%ylim([-0.1 0.3]);
vline(0,'b'); vline(laseron/10,'b');
box off; hold off;
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'FontSize',8,'FontName','Arial');
print(align_fig,[image_dest 'dfOvF_btmbase_align_laseron.pdf'],'-r600','-dpdf');
savefig([image_dest 'dfOvF_btmbase_align_laseron.fig']);
% print(align_fig,[image_dest 'dfOvFavglaseroff_align_laseron_frames.pdf'],'-r600','-dpdf');
% savefig([image_dest 'dfOvFavglaseroff_align_laseron_frames.fig']);
