%paths
ds = 'YM90K_DART_exptList';
eval(ds)
iexp = 7;
for day = 1:2

    mouse = expt(iexp).mouse;
    date = expt(iexp).dates{day};
    run = expt(iexp).runs{day};
    time = expt(iexp).times{day};
    
    data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\ACh';
    %CD = [data_pn '\Data\2P_images\' mouse '\' date '\' run];
    CD = [data_pn '\Data\2P_data\' mouse '\' date '\' run];
    cd(CD);
    fn = [run '_000_000_eye.mat'];
    
    %load data
    data_temp = load(fn);
    data_temp = squeeze(data_temp.data);
    
    %crop frames to match mworks data
    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' mouse '-' date '-' time '.mat'];
    load(fName);
    nFrames = input.counterValues{end}(end);
    nTrials = size(input.counterValues,2);
    data = data_temp(:,:,1:nFrames);      % the raw images...
    
    % Crop image to isolate pupil 
    %(bright spots can be mistaken for pupil)
    if day == 1
        [data_crop rect] = cropEyeData(data);
    else
        [data_crop rect] = cropEyeData(data,rect);
    end
    
    % measure pupil position/diameter
    rad_range = [3 15]; %adjust to expected range of pupil size (if low end is too small then may find noisy bright stuff)
    Eye_data = extractEyeData(data_crop,rad_range);
    %if pupil not found reliably, adjust the image cropping or the rad_range
    
    % align to stimulus presentation
    [rad centroid] = alignEyeData(Eye_data,input);
    
    % wheel data
    wheel_data = wheelSpeedCalc(input,32,'purple');
    wheel_trial = mean(reshape(wheel_data,[input.nScansOff+input.nScansOn nTrials]),1);
    data_out = fullfile(data_pn,'Analysis','2p_analysis',mouse,date,run);
    save(fullfile(data_out,'pupil.mat'),'rect','rad_range','Eye_data','rad','centroid','wheel_trial')

end
%% 

data_pn = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\ACh';
figure;
for id = 1:2
    data_out = fullfile(data_pn,'Analysis','2p_analysis',mouse,expt(iexp).dates{id},expt(iexp).runs{id});
    load(fullfile(data_out,'pupil.mat'))
    subplot(2,2,1)
    scatter(wheel_trial,rad.stim)
    hold on
    subplot(2,2,2)
    histogram(rad.stim,[0:0.02:0.5])
    hold on
    subplot(2,2,3)
    histogram(rad.stim(find(wheel_trial>2)),[0:0.02:0.5])
    hold on
    subplot(2,2,4)
    histogram(rad.stim(find(wheel_trial<2)),[0:0.02:0.5])
    hold on
end
subplot(2,2,1)
xlabel('Running speed')
ylabel('Pupil diameter')
legend(expt(iexp).dates)
subplot(2,2,2)
xlabel('Pupil diameter')
ylabel('Count')
subplot(2,2,3)
xlabel('Pupil diameter')
ylabel('Count')
title('Running')
subplot(2,2,4)
xlabel('Pupil diameter')
ylabel('Count')
title('Stationary')
suptitle(mouse)
data_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\ForCeline\pupil';
print(fullfile(data_out,[mouse '_' expt(iexp).dates{1} '-' expt(iexp).dates{2} '_dirXcon_pupil.pdf']),'-dpdf','-bestfit')

figure;
for id = 1:2
    data_out = fullfile(data_pn,'Analysis','2p_analysis',mouse,expt(iexp).dates{id},expt(iexp).runs{id});
    load(fullfile(data_out,'pupil.mat'))
    run_trial = find(wheel_trial>2);
    ntrial = length(wheel_trial);
    subplot(2,1,id)
    scatter(1:ntrial,rad.stim)
    hold on
    scatter(run_trial,rad.stim(run_trial))
    title(expt(iexp).dates{id})
    xlabel('Trial')
    ylabel('Pupil diameter')
    legend('stationary','running','location','southwest')
    ylim([0 0.6])
end
suptitle(mouse)
data_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\ForCeline\pupil';
print(fullfile(data_out,[mouse '_' expt(iexp).dates{1} '-' expt(iexp).dates{2} '_dirXcon_pupilTC.pdf']),'-dpdf','-bestfit')
    
            
