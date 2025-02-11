% =============================================================
% Simple helper functions and MATLAB structures demonstrating
% how to read and manipulate SpikeGLX meta and binary files.
%
% The most important part of the demo is ReadMeta().
% Please read the comments for that function. Use of
% the 'meta' structure will make your data handling
% much easier!
%
function [Waveforms, SampleTS, reporter] = LongTraceReadSGLXData(dataLength, n, timestamps, chan)
% dataLength = how much data to extract (waveform length)
% n = how many waveforms to extract
% timestamps = when unit is firing
% ch = channel in SPIKEGLX NOTATION- zero-based

% Ask user for binary file
[binName,path] = uigetfile('*.bin', 'Select Binary File');

% Parse the corresponding metafile
meta = ReadMeta(binName, path);

chan = chan +1;
step = floor(length(timestamps)/n);

timestamps = ((timestamps)); % catch the begining of the waveform
nSamp = floor(dataLength * SampRate(meta));
nSamp = int64(nSamp);
Waveforms = zeros(nSamp, n);
SampleTS = zeros(n,1);


for i = 1:n
  randn = randi([(((i-1)*step)+1),(i*step)]);
SampleTS(i) = timestamps(randn); 
samp0 = SampleTS(i)*SampRate(meta);    
samp0 = int64(samp0);
dataArray = ReadBin(samp0, nSamp, meta, binName, path);

% % % change y values to volts
% For an analog channel: gain correct saved channel ch (1-based for MATLAB).
ch = chan;

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;

    if strcmp(meta.typeThis, 'imec')
        dataArray = GainCorrectIM(dataArray, [ch], meta);
    else
        dataArray = GainCorrectNI(dataArray, [ch], meta);
    end
  % % % 

Waveforms(:,i) = dataArray(chan,:);
t = 1:1:((length(Waveforms(:,1))));
t = t/SampRate(meta) + SampleTS(i);
t = t.';
size(t)
size(Waveforms(:,i))
end   
figure
% % % plot(t, Waveforms(:,i));
plot(Waveforms(:,i));
%xlim([t(1) t(end)]);
hold on

TSwindow = (timestamps(( timestamps >= SampleTS(i)) & (timestamps <=(SampleTS(i) + dataLength))));
mininSamp = floor(.001*SampRate(meta));
mininSamp = int64(mininSamp)
miniWaveforms = zeros(mininSamp, length(TSwindow));
miniTS = 1:1:floor(.001*int64(SampRate(meta)));
miniTS = miniTS.';

    for i = 1:(length(TSwindow)-1)
    miniTS = [miniTS miniTS];
    end

%reporter = SampleTS(i);

for i = 1:length(TSwindow)
    %figure
samp0 = int64(TSwindow(i)*SampRate(meta));
dataArray1 = ReadBin(samp0, mininSamp, meta, binName, path);

% % % change y values to volts
% For an analog channel: gain correct saved channel ch (1-based for MATLAB).
ch = chan;

% For a digital channel: read this digital word dw in the saved file
% (1-based). For imec data there is never more than one saved digital word.
dw = 1;

    if strcmp(meta.typeThis, 'imec')
        dataArray1 = GainCorrectIM(dataArray1, [ch], meta);
    else
        dataArray1 = GainCorrectNI(dataArray1, [ch], meta);
    end
  % % % 
miniWaveforms(:,i) = dataArray1(chan,:);
%reporter = miniWaveforms;
%TSwindowI = (TSwindow(i));

%miniTS = ((TSwindowI - TSwindow(1))*int64(SampRate(meta))) : 1 : floor(.001*int64(SampRate(meta))-1); 
 %?
%reporter = miniTS;
miniTS(:,i) = miniTS(:,i) + (TSwindow(i)- TSwindow(1))*SampRate(meta);
reporter = miniTS;
% % %t = (double(miniTS(:,i))/(SampRate(meta))) + (TSwindow(i));
%reporter = miniTS;
%reporter = TSwindowI;
%plot(miniTS, miniWaveforms(:,i), 'k');
%miniTS = miniTS + TSwindow(i)+int64(SampRate(meta))
plot(miniTS(:,i), miniWaveforms(:,i), 'k');
% % %plot(t, miniWaveforms(:,i), 'k');
%figure
%plot(miniWaveforms(:,i), 'k');
%reporter = dataArray1;
end



end % DemoReadSGLXData


% =========================
% General Utility Functions
% =========================


% =========================================================
% Parse ini file returning a structure whose field names
% are the metadata left-hand-side tags, and whose right-
% hand-side values are MATLAB strings. We remove any
% leading '~' characters from tags because MATLAB uses
% '~' as an operator.
%
% If you're unfamiliar with structures, the benefit
% is that after calling this function you can refer
% to metafile items by name. For example:
%
%   meta.fileCreateTime  // file create date and time
%   meta.nSavedChans     // channels per timepoint
%
% All of the values are MATLAB strings, but you can
% obtain a numeric value using str2double(meta.nSavedChans).
% More complicated parsing of values is demonstrated in the
% utility functions below.
%
function [meta] = ReadMeta(binName, path)

    % Create the matching metafile name
    [dumPath,name,dumExt] = fileparts(binName);
    metaName = strcat(name, '.meta');

    % Parse ini file into cell entries C{1}{i} = C{2}{i}
    fid = fopen(fullfile(path, metaName), 'r');
% -------------------------------------------------------------
%    Need 'BufSize' adjustment for MATLAB earlier than 2014
%    C = textscan(fid, '%[^=] = %[^\r\n]', 'BufSize', 32768);
    C = textscan(fid, '%[^=] = %[^\r\n]');
% -------------------------------------------------------------
    fclose(fid);

    % New empty struct
    meta = struct();

    % Convert each cell entry into a struct entry
    for i = 1:length(C{1})
        tag = C{1}{i};
        if tag(1) == '~'
            % remake tag excluding first character
            tag = sprintf('%s', tag(2:end));
        end
        meta = setfield(meta, tag, C{2}{i});
    end
end % ReadMeta


% =========================================================
% Read nSamp timepoints from the binary file, starting
% at timepoint offset samp0. The returned array has
% dimensions [nChan,nSamp]. Note that nSamp returned
% is the lesser of: {nSamp, timepoints available}.
%
% IMPORTANT: samp0 and nSamp must be integers.
%
function dataArray = ReadBin(samp0, nSamp, meta, binName, path)

    nChan = str2double(meta.nSavedChans);

    nFileSamp = str2double(meta.fileSizeBytes) / (2 * nChan);
    samp0 = max(samp0, 0);
    %samp0
    %nSamp = samp0+nSamp;
    nSamp = min(nSamp, nFileSamp - samp0);
    %nSamp
    %nSamp-samp0 

    sizeA = [nChan, nSamp];

    fid = fopen(fullfile(path, binName), 'rb');
    fseek(fid, samp0 * 2 * nChan, 'bof');
    dataArray = fread(fid, sizeA, 'int16=>double');
    fclose(fid);
end % ReadBin




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
% Return a multiplicative factor for converting 16-bit
% file data to voltage. This does not take gain into
% account. The full conversion with gain is:
%
%   dataVolts = dataInt * fI2V / gain.
%
% Note that each channel may have its own gain.
%
function fI2V = Int2Volts(meta)
    if strcmp(meta.typeThis, 'imec')
        fI2V = str2double(meta.imAiRangeMax) / 512;
    else
        fI2V = str2double(meta.niAiRangeMax) / 32768;
    end
end % Int2Volts


% =========================================================
% Return array of original channel IDs. As an example,
% suppose we want the imec gain for the ith channel stored
% in the binary data. A gain array can be obtained using
% ChanGainsIM() but we need an original channel index to
% do the look-up. Because you can selectively save channels
% the ith channel in the file isn't necessarily the ith
% acquired channel, so use this function to convert from
% ith stored to original index.
%
% Note: In SpikeGLX channels are 0-based, but MATLAB uses
% 1-based indexing, so we add 1 to the original IDs here.
%
function chans = OriginalChans(meta)
    if strcmp(meta.snsSaveChanSubset, 'all')
        chans = (1:str2double(meta.nSavedChans));
    else
        chans = str2num(meta.snsSaveChanSubset);
        chans = chans + 1;
    end
end % OriginalChans


% =========================================================
% Return counts of each imec channel type that compose
% the timepoints stored in binary file.
%
function [AP,LF,SY] = ChannelCountsIM(meta)
    M = str2num(meta.snsApLfSy);
    AP = M(1);
    LF = M(2);
    SY = M(3);
end % ChannelCountsIM

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


% =========================================================
% Return gain for ith channel stored in the nidq file.
%
% ichan is a saved channel index, rather than an original
% (acquired) index.
%
function gain = ChanGainNI(ichan, savedMN, savedMA, meta)
    if ichan <= savedMN
        gain = str2double(meta.niMNGain);
    elseif ichan <= savedMN + savedMA
        gain = str2double(meta.niMAGain);
    else
        gain = 1;
    end
end % ChanGainNI


% =========================================================
% Return gain arrays for imec channels.
%
% Index into these with original (acquired) channel IDs.
%
function [APgain,LFgain] = ChanGainsIM(meta)

    if isfield(meta,'imDatPrb_dock')
        [AP,LF,~] = ChannelCountsIM(meta);
        % NP 2.0; APgain = 80 for all channels
        APgain = zeros(AP,1,'double');
        APgain = APgain + 80;
        % No LF channels, set gain = 0
        LFgain = zeros(LF,1,'double');
    else
        % 3A or 3B data?
        % 3A metadata has field "typeEnabled" which was replaced
        % with "typeImEnabled" and "typeNiEnabled" in 3B.
        % The 3B imro table has an additional field for the
        % high pass filter enabled/disabled
        if isfield(meta,'typeEnabled')
            % 3A data
            C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d', ...
                'EndOfLine', ')', 'HeaderLines', 1 );
        else
            % 3B data
            C = textscan(meta.imroTbl, '(%*s %*s %*s %d %d %*s', ...
                'EndOfLine', ')', 'HeaderLines', 1 );
        end
        APgain = double(cell2mat(C(1)));
        LFgain = double(cell2mat(C(2)));
    end
end % ChanGainsIM


% =========================================================
% Having acquired a block of raw nidq data using ReadBin(),
% convert values to gain-corrected voltages. The conversion
% is only applied to the saved-channel indices in chanList.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [1:MN]      % all MN chans (MN from ChannelCountsNI)
%   [2,6,20]    % just these three channels
%
function dataArray = GainCorrectNI(dataArray, chanList, meta)

    [MN,MA] = ChannelCountsNI(meta);
    fI2V = Int2Volts(meta);

    for i = 1:length(chanList)
        j = chanList(i);    % index into timepoint
        conv = fI2V / ChanGainNI(j, MN, MA, meta);
        dataArray(j,:) = dataArray(j,:) * conv;
    end
end


% =========================================================
% Having acquired a block of raw imec data using ReadBin(),
% convert values to gain-corrected voltages. The conversion
% is only applied to the saved-channel indices in chanList.
% Remember saved-channel indices are in range [1:nSavedChans].
% The dimensions of the dataArray remain unchanged. ChanList
% examples:
%
%   [1:AP]      % all AP chans (AP from ChannelCountsIM)
%   [2,6,20]    % just these three channels
%
function dataArray = GainCorrectIM(dataArray, chanList, meta)

    % Look up gain with acquired channel ID
    chans = OriginalChans(meta);
    [APgain,LFgain] = ChanGainsIM(meta);
    nAP = length(APgain);
    nNu = nAP * 2;

    % Common conversion factor
    fI2V = Int2Volts(meta);

    for i = 1:length(chanList)
        j = chanList(i);    % index into timepoint
        k = chans(j);       % acquisition index
        if k <= nAP
            conv = fI2V / APgain(k);
        elseif k <= nNu
            conv = fI2V / LFgain(k - nAP);
        else
            continue;
        end
        dataArray(j,:) = dataArray(j,:) * conv;
    end
end


