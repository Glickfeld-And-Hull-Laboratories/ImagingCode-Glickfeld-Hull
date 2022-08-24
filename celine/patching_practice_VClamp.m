%% looking at voltage clamp data
wd = 'Z:\home\celine\Data\patch_data';
cd(wd)
%add something to change to the directory for thie data I currently want,
%by date
[Idata,~,Iparams] = abfload('22808000.abf');
%samples X one channel (data) X sweeps


nSweeps = size(Idata,3);
sf=20000; %enter sampling rate


%% identifies the stimulus time, duration, and magnitude for each sweep
onThreshold = 5;
eStim = (squeeze(Idata(:,2,:)));
eStim = (eStim-min(eStim))/10;
temp_stimOn =eStim>onThreshold;
stimParams=struct;
for iSweep = 1:nSweeps
    stimParams.stimOn(iSweep) = find(temp_stimOn(:,iSweep),1,'first');
    stimParams.stimOff(iSweep)=(find(temp_stimOn(:,iSweep),1,'last')+1);
    stimParams.stimMag(iSweep)=max(eStim(:,iSweep))-min(eStim(:,iSweep));
    bsln=median(Idata((stimParams.stimOn(iSweep)-100):(stimParams.stimOn(iSweep)-50),1,iSweep));
    peak=max(Idata((stimParams.stimOn(iSweep)-50):(stimParams.stimOn(iSweep)+50),1,iSweep));
    stimParams.R_access(iSweep) = (stimParams.stimMag(iSweep)/(peak-bsln))*1000;
    clear bsln peak
    stimParams.stimDur = stimParams.stimOff(iSweep)-stimParams.stimOn(iSweep);
end


% the access resistance
figure;plot(stimParams.R_access);xlabel("sweep number");ylabel("R access")%to see how access changes over the sweeps
%% Plots everything 
t = make_time(Idata(:,:,:),sf,1);

figure;
subplot(2,1,1)
hold on
for iSweep = 1:nSweeps
    plot(t, Idata(:,1,iSweep));
end


%find the time of the test pulse
%find the access resistance

subplot(2,1,2)
hold on
for iSweep = 1:nSweeps
    plot(t, eStim(:,iSweep));
end
hold off 
%% 






