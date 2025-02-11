% Uses txt files extracted with CatGT (XA=0,.5,.5,0 and iXA=0,.5,.5,0: Look carefully in upper left corner of dialog box) to
% create states for timestamps for running data. Can be put through TPrime
% but probably not nessesary because I remove one second on either side of
% trial change timestamps.
%
% plots a matrix (time) with time (sec) in the first column and corresponding
% states (1= still, 2 = running) in second column.
%
% creates ifrun, a 2-column matrix where the first column is the start time
% of every trial, starting at 0. The second column has a 1 if the trial is
% stationary and a 2 if the trial is locomotion. Written to analyze
% Shuyangs running data. Works with FRstructRunNorun to analyze Shuyangs
% running data
%
% MEH 6/28/21

function ifrun = ifRun();

[txtname,path] = uigetfile('*.txt', 'Select XAfile ');
fid = fopen(txtname);
UpTTL = fscanf(fid, '%f');
UpTTL(:,2) = 2;
[txtname,path] = uigetfile('*.txt', 'Select iXAfile ');
fid = fopen(txtname);
DownTTL = fscanf(fid, '%f');
DownTTL(:,2) = 1;


TTL = [UpTTL; DownTTL];
sortTTL = sortrows(TTL);
firstRow = [0 1];
TTL = [firstRow; sortTTL];
plot(TTL(:,1), TTL(:,2));

time = (0:1:TTL(end,1)+30).';

for k = 1:(length(TTL))
    state = TTL(k,2);
    index = time(:,1)<TTL(k+1);
    lasttimeindex = length(find(index));
    index= time(:,1)<TTL(k);
    firsttimeindex = length(find(index));
    if firsttimeindex == 0
        firsttimeindex =1;
    end
    time(firsttimeindex:lasttimeindex,2) = state;
end

%put in state for the last 30 seconds
index = find(time(:,2) == 0);
firsttimeindex = index(1);
lasttimeindex = index(end);
time(firsttimeindex:lasttimeindex,2) = state;

% make bar-chart like graph showing running and not running times
area(time(:,1), time(:,2));
ylim([1 3]);

% create ifrun, a 2-column matrix where the first column is the start time
% of every trial, starting at 0. The second column has a 1 if the trial is
% stationary and a 2 if the trial is locomotion. Written to analyze
% Shuyangs running data. Works with FRstructRunNorun to analyze Shuyangs
% running data

%find start time, for trialLength long trials (seconds)
trialLength = 30;
maybestart = TTL(2,1);

while maybestart > 30
    maybestart = maybestart - 30;
end
start = maybestart;

trialStart = start;
m = 1;
lastTrialStart = time(end,1);
while (trialStart < lastTrialStart)
ifrun(m,1) = trialStart;
checktimeIndex = floor(trialStart + 5);
ifrun(m,2) = time((checktimeIndex), 2);
m = m+1;
trialStart = trialStart + trialLength;
end