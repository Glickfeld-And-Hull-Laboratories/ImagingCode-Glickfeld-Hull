function WaveformsAdj = reAlignWF(Waveforms, SampleTS, time, preTime, SampRate, peakWin, n);

AvgWvF = avgeWaveforms(Waveforms);
prom = (max(AvgWvF)-min(AvgWvF))/4;

peakWin;
winMin = round(peakWin(1));
winMax = round(peakWin(2));

    for p = 1:n
        ThisWave = Waveforms(:,p);
        FindShift = ThisWave;
    %plot(time, FindShift, 'k');
    FindShift(1:winMin)=NaN;
    FindShift(winMax:end)=NaN;
    TF = islocalmin(FindShift, 'MinProminence', prom);
    %TF = islocalmin(FindShift);
    k = find(TF); %the index for the minima associated with the negative spike
    TimeZeroIn = preTime * SampRate; % the index where we want to align all the minima
    TimeZeroIn = double(TimeZeroIn);
    shift = k - TimeZeroIn
     if length(k) > 1 % if there are two minima in the window outstide of the Minimum Prominience, take the one that is closer to the timestamp of the spike 
       [m,I] = min(abs(shift));
       shift = shift(I)
    end
    if shift > 0
        ThisWave = ThisWave((shift+1):end);
        ThisWave = [ThisWave; NaN(shift, 1)];
    end
    if shift < 0
        ThisWave = ThisWave(1:(end+shift));
        ThisWave = [NaN(-shift, 1); ThisWave];
    end
    WaveformsAdj(:, p) = ThisWave;
    %reporter = TF;
    %TF = islocalmin(ThisWave, 'MinProminence', prom);
    %plot(time, ThisWave,time(TF),ThisWave(TF),'r*');
    end
   % TF = islocalmin(Waveforms(:,i));
    %plot(x,A,x(TF),A(TF),'r*')

AvgWvF = avgeWaveforms(Waveforms);
prom = (max(AvgWvF)-min(AvgWvF))/2;
TF = islocalmin(AvgWvF, 'MinProminence', prom);
%figure
%plot(time,AvgWvF,time(TF),AvgWvF(TF),'r*')
end