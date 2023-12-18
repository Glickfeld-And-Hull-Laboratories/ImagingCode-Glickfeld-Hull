function [baselineCounts, drugCounts, baselineAmps, drugAmps] = eventCountAndAmpDRUG(dataIN,drugStart,totalSweeps)
% Function that takes a file name for a .csv file of events and the sweep at which the drug arrived, and
% returns four values: the mean number of events per sweep in the baseline
% condition and in the drug condition, and the mean event amplitude in each
% condition. 

    % I'm allowing a maximum amplitude of 0.09 nAmp and a maximum rise
    % slope of 0.05 nA/mS. Events above these are values are removed.
    ampThreshold = 0.09;
    slopeThreshold = 0.1;

    drugBuffer=drugStart+10; %10 sweeps to allow drug to permeate slice
    drugDuration = 20; %20 sweeps assessed with the drug on
    %determine the last sweep based on how many sweeps I have - it will 
    % either be after 20 sweeps with the drug on or at the end of the
    % recording, whichever is less
    endSweep=min([totalSweeps,(drugBuffer+drugDuration)]);

    sweepLength = 9.5; %this is the sweep duration I usually use
    
    path = fullfile('\home','celine','Analysis','ePhys');
    cd(path) %go to where the data files are
    inputFile = [dataIN,'.csv'];
    fullData =readtable(inputFile); %open the desired file as a table
    %find the column number for the peak amplitude
    traceCol = find(string(fullData.Properties.VariableNames) == "Trace");
    peakCol = find(string(fullData.Properties.VariableNames) == "Peak");
    riseSlopeCol = find(string(fullData.Properties.VariableNames) == "MaxRiseSlope_nAperMs");
    %take the columns I need and convert to an array
    data=fullData(:,[traceCol peakCol riseSlopeCol]);
    data=table2array(data);
    %now 1=Trace, 2=peak, and 3=rise slope
    %eliminate events the fail cirteria

    %data(find(abs(data(:,2))>0.07),:)=[];
    mean(abs(data(:,2)));
    ampCutOFf = abs(data(:,2))>ampThreshold;
    data(ampCutOFf,:)=[];
    mean(abs(data(:,2)));
    
    mean(abs(data(:,3)));
    slopeCutOFf = abs(data(:,3))>slopeThreshold;
    data(slopeCutOFf,:)=[];
    mean(abs(data(:,3)));
    %endSweep, the total number of sweeps to interrogate, will be set to the
    %end of the drug period that we're interested in
   
    counts=zeros(endSweep,1);
    amps=nan(endSweep,1);
    
    for iSweep = 1:endSweep
        if length(find(data(:,1)==iSweep))>0
            amps(iSweep)=abs(nanmean(data(find(data(:,1)==iSweep),2)));
            counts(iSweep)=length(find(data==iSweep));
        end
    end
    
    counts = (counts./sweepLength)*60;
    baselineCounts = nanmean(counts(1:drugStart));
    drugCounts = nanmean(counts(drugBuffer:endSweep));
    baselineAmps = nanmean(amps(1:drugStart));
    drugAmps = nanmean(amps(drugBuffer:endSweep));


end