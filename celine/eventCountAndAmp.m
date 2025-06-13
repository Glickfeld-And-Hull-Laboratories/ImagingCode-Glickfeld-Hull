function [meanCounts, meanAmps,nTotal, nRejectedHigh,nRejectedLow] = eventCountAndAmp(dataIN,startSweep,stopSweep,conversionFact)
% Function that takes a file name for a .csv file of events and the sweep at which the drug arrived, and
% returns four values: the mean number of events per min and the mean
% amplitude

    % I'm allowing a maximum amplitude of 60 pAmp and a maximum rise
    % slope of 0.1 nA/mS. Events above these are values are removed.
    ampMax = 175;
    ampMin = 15;
    slopeThreshold = 1000;

    sweepLength = 9.5; %this is the sweep duration I usually use
    
    path = fullfile('\home','celine','Analysis','ePhys');
    cd(path) %go to where the data files are
    inputFile = [dataIN,'.csv'];
    fullData =readtable(inputFile); %open the desired file as a table
    %take the columns I need and convert to an array
    data=fullData(:,[1 8 20]);
    data=table2array(data);
    data(:,2:3)=data(:,2:3)*conversionFact;
    %now 1=Trace, 2=peak, and 3=rise slope
    %eliminate events the fail cirteria

   
    ampTooHigh = abs(data(:,2))>ampMax;
    nRejectedHigh = sum(ampTooHigh);
    data(ampTooHigh,:)=[];

    ampTooLow = abs(data(:,2))<ampMin;
    nRejectedLow = sum(ampTooLow);
    data(ampTooLow,:)=[];
    
    nTotal = size(data,1);
    
    slopeCutOFf = abs(data(:,3))>slopeThreshold;
    data(slopeCutOFf,:)=[];


    %endSweep, the total number of sweeps to interrogate, will be set to the
    %end of the drug period that we're interested in
    nSweeps=stopSweep-startSweep+1;
    counts=zeros(nSweeps,1);
    amps=nan(nSweeps,1);
    
    for iSweep = startSweep:stopSweep
        if length(find(data(:,1)==iSweep))>0
            amps(iSweep)=abs(nanmean(data(find(data(:,1)==iSweep),2)));
            counts(iSweep)=length(find(data==iSweep));
        end
    end
    
    counts = (counts./sweepLength)*60;
    meanCounts = nanmean(counts);
    if meanCounts ==0
        meanAmps=NaN;
    else
    meanAmps = nanmean(amps);
    end


end