%% ephys practice

%load ABF file
% --> generalize using file list
[Vdata,~,Vparams] = abfload('22718000.abf');
%samples X channels (1=data, 2=stim) X sweeps/trials

sf=20000; %enter sampling rate

nSweeps = size(Vdata,3);
onThreshold = 0.3;
eStim = abs(squeeze(Vdata(:,2,:)));
temp_stimOn =(eStim./max(eStim))>onThreshold;
for iSweep = 1:nSweeps
    stimOn(iSweep) = find(temp_stimOn(:,iSweep),1,'first');
    stimOff(iSweep)=(find(temp_stimOn(:,iSweep),1,'last')+1);
end

stimDur = stimOff-stimOn;

%currently this doesn't account for spiking
%input resistance seems to change a decent amount depending on current step

%make a matrix to save these values
for iSweep = 1:nSweeps
    Ibsln(iSweep)=mean(Vdata(1:stimOn,2,iSweep)); %most likely to be approx 0, but check
    Istep(iSweep)=mean(Vdata(stimOn:stimOff,2,iSweep))-Ibsln(iSweep);
    Vbsln(iSweep) =mean(Vdata(1:stimOn,1,iSweep)); %resting membrane potential
    Vstep(iSweep)=abs(mean(Vdata(stimOn:stimOff,1,iSweep))-Vbsln(iSweep));
    %use V=IR to find R, Vstep/Istep = R, and multiply by 1000 to convert to
    %megaOhms

end




%%

%to find the input restance based on average of the -50 pA current step
%steps
temp_step = find(abs(Istep-(-50))==min(abs(Istep-(-50)))); %find the current step closest to -/+50
Rinput=abs((Vstep(temp_step))/abs(Istep(temp_step)))*1000;

clear temp_step
scatter(Istep,Vstep)
%%

% Plots everything 
t = make_time(Vdata(:,:,:),sf,1);

figure;
subplot(2,1,1)
hold on
for iSweep = 8:8
    plot(t, Vdata(:,1,iSweep));
end
hold off 

subplot(2,1,2)
hold on
for iSweep = 8:8
    plot(t, Vdata(:,2,iSweep));
end
hold off 
%% find peak voltage

peakAmps = []; % 
for iSweep = 1:nSweeps
    maxAmp = max(Vdata(stimOn:stimOff,1,iSweep))
    peakAmps = [peakAmps;maxAmp]
end 



%% Jenny's clean data function + my function to find spike times and rate
%I don't understand what the points are here
notZero=logical(stimOn~=1); %identify sweeps with no current step
Vclean = cell(1,nSweeps);
FR_by_step = zeros(2,nSweeps);


%need to skip the 0 step
for iSweep = 1:nSweeps
    if notZero(iSweep)
        temp_data=cleanDataCeline(Vdata(:,:,iSweep),1,stimOn(iSweep)-50); 
        Vclean{iSweep}=temp_data;
        [FR_by_step(2,iSweep),Vclean{iSweep}.spikeTimes] = SpikeRate(Vclean{iSweep}.pts,sf);
        FR_by_step(1,iSweep) = Istep(iSweep);
    end
end

%plot the relationship of firing rate and current step size
figure;scatter(FR_by_step(1,:),FR_by_step(2,:));xlabel("Current step");ylabel("Firing rate");




