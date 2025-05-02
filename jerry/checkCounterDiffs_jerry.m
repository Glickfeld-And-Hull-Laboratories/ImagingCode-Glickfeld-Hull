clear all; clear global;  
clc

ds = 'DART_expt_info'; %dataset info
dataStructLabels = {'contrastxori'};
rc = behavConstsDART; %directories
eval(ds);
doGreenOnly = true;
doCorrImg = true;

day_id = input('Enter day id ');

if computer == 'GLNXA64'
    beh_prefix = strcat(isilonName,'/All_Staff/Behavior/Data/');
elseif string(hostname) == 'NB-NUKE'
    beh_prefix = strcat('Z:/All_Staff/Behavior/Data/');
else
    beh_prefix = strcat('Z:\Behavior\Data\');
end

mouse = expt(day_id).mouse;
expDate = expt(day_id).date;
dat = 'data-';
time = eval(['expt(day_id).' cell2mat(dataStructLabels) '_time']);
tHostname = lower(hostname);
    [s,tUsername] = dos('ECHO %USERNAME%');
    switch tHostname
        case {'nb-nuke'}
            if username == 'cc735' 
                fName = ['Z:\All_staff\Behavior\Data\' dat mouse '-' expDate '-' time{1} '.mat'];

            else
                fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' dat mouse '-' expDate '-' time{1} '.mat'];
            end
        case{'nb-hubel'}
                if username == 'cc735'
                    fName = ['Z:\Behavior\Data\' dat mouse '-' expDate '-' time{1} '.mat'];
                else
                    fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\' dat mouse '-' expDate '-' time{1} '.mat'];
                 end
    end
    
load(fName); %load the mworks behavioral file

plotCellArrayDifferencesTwoInputs(input.counterTimesUs,input.counterValues)
sgtitle(['Mouse ' mouse ', ' expDate])