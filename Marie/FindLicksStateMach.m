% =============================================================
% Adapted from DemoReadSGLXData that came from SpikeGLX materials from Bill
% Karsh. Takes inputs of Juice delivery timestamps (extracted using CatGT)
% and and exported nidaq channel that has only analog4 (licking channel)
% and all digiital channels (could proably use only analog4 now, and would
% like to modify to use original nidaq channel). Export is made in SpikeGLX
% file viewer.
%
% FindLicks uses tmin and tmax to create a window around each juice
% delivery event and extract licks before (tmin) and after (tmax) juice
% delivery.
%
% FindLicks returns a cell array with one field containing each juice time
% and the other containing a vector with associated licktimes (tmin before
% and tmax after each delivery). 3rd field minsample 4th field maxsample
% 5th field dataArray.
%
% MEH 2/25/21
%


function [JuiceLicks] = FindLicksStateMach(JuiceTimes, tmin, tmax, level1, level2)
%clear
%threshold = 1; % change this to change threshold for detecting licks
% Ask user for binary file
[binName,path] = uigetfile('*.bin', 'Select Binary File');
fprintf ('File chosen is %s ', path);
fprintf ('%s \n', binName);
tic 

% Parse the corresponding metafile
meta = ReadMeta(binName, path);

fTime = str2double(meta.fileTimeSecs);

JuiceLicks = cell(length(JuiceTimes),5);
        
for m = 1:length(JuiceTimes) % change to length(JuiceTimes) after testing
    
    fprintf('\n JuiceTime = %d \n', round(JuiceTimes(m)));
    fprintf('JuiceIndex = %d \n', m);
    toc;
    
    readMin = int64((JuiceTimes(m)- tmin )*SampRate(meta));
    if readMin < 0
        readMin = 0;
    end
                        %JuiceLicks{m,3}= readMin; % for troubleshooting

    readMax = int64((JuiceTimes(m) + tmax )*SampRate(meta));
    if readMax > fTime*SampRate(meta)
         readMax = fTime*SampRate(meta);
    end
                        %JuiceLicks{m,4}=readMax; % for troubleshooting

    dataArray = ReadBin(readMin, readMax, meta, binName, path);

    % For an analog channel: gain correct saved channel ch (1-based for MATLAB).
        ch = 1;

    dataArray = GainCorrectNI(dataArray, [ch], meta);
    
    dataArray = dataArray(1,:); %If the exported data has different number or configuration of channels, this will need to be adjusted


    JuiceLicks{m,3} = dataArray;
    % [smArray, ~] = smoothdata(JuiceLicks{m,3},'gaussian', 70000); %get dataArray from MEHreadGLX and smooth it

    thresholdData = dataArray > level1;
            %JuiceLicks{m,5} = thresholdData; % for troubleshooing
          
    indexOnes = find(thresholdData);
    licks = zeros(1,length(thresholdData));
    state = 'waiting';
    for i = 2:(length(indexOnes))                       % for every timepoint where the trace is above the threshold, check if it is the first trace above the threshold.
         % for every timepoint where the trace is above the threshold, check if it is the first trace above the threshold.
            if (~(thresholdData((indexOnes(i)-1))))     % if it's just crossing at i
                sample = indexOnes(i);
                down = find(~thresholdData((sample+1):end),1);     % find out where it drops back down
                stop = (down + sample);                   % set sample where it drops below thresh1 as stop
                win = [sample; stop];                        % window covers data from coming above threshold until it goes below again
                
                if (stop - sample) > SampRate(meta)*.002     % if it stays above threshold for only a short time, go back to waiting. If it stays above thresh1, lets check for thresh2 in window
                    state = 'thresh1';
                else                                    
                    state = 'waiting2';
                end
                if  strcmp(state, 'thresh1')            % If we are above thresh1 and it's potentially a real hit
                    thresholdData2 = dataArray > level2;    %find out if trace goes aboce thresh2 during the window.
                    if any(thresholdData2(sample:stop))             %so, validate point i as begining a lick
                        licks(1,sample) = 1;                       % If so, mark that timepoint with a 1 in output data
                        state = 'thresh2';
                    else
                        state = 'waiting';               % If you get aboce thresh1 but not thresh2 during the window, go back to waiting
                    end
                    state = 'waiting';                   % once you count the lick, go back to waiting
                end
            end
    end
                        %JuiceLicks{m,7} = licks;
    clear LickTimes;
    LickTimes = (find(licks)/SampRate(meta))+ JuiceTimes(m)-tmin;
         if (~isempty(LickTimes))
             JuiceLicks{m,2} = LickTimes;
         else
             LickTimes = [NaN];
               fprintf('no licks for %d \n', round(JuiceTimes(m)));
               JuiceLicks{m,2} = LickTimes;
         end

JuiceLicks{m,1} = JuiceTimes(m);

end
%licktimes = LickTimes;
%assignin('base', 'LickTimes', LickTimes);
end % DemoReadSGLXData

% ====================
% 

% =========================
% General Utility Functions
% =========================











% =========================================================
% Return sample rate as double.
%
function srate = SampRate(meta)
    if strcmp(meta.typeThis, 'imec')
        srate = str2double(meta.imSampRate);
    else
        srate = str2double(meta.niSampRate);
    end
end % SampRate

% =========================================================
% Return file time as double.
%
function fTime = FileTime(meta)
    if strcmp(meta.typeThis, 'imec')
        fTime = str2double(meta.fileTimeSecs);
    else
        fTime = str2double(meta.fileTimeSecs);
    end
end % fTime






% =========================================================
% Return counts of each nidq channel type that compose
% the timepoints stored in binary file.
%
function [MN,MA,XA,DW] = ChannelCountsNI(meta)
    M = str2num(meta.snsMnMaXaDw);
    MN = M(1);
    MA = M(2);
    XA = M(3);
    DW = M(4);
end % ChannelCountsNI


