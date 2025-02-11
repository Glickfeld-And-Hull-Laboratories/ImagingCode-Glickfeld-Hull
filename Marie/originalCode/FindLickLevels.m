function FindLickLevels(JuiceTimes, tmin, tmax, n)
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

Prelim = cell(length(JuiceTimes),5);
for m = 1:n  % replace 'n' with 'length(JuiceTimes' to go back to function FindLicksStateMach % change to length(JuiceTimes) after testing
    
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


    Prelim{m,3} = dataArray;
    % [smArray, ~] = smoothdata(JuiceLicks{m,3},'gaussian', 70000); %get dataArray from MEHreadGLX and smooth it
    plot(dataArray);

  
 



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


