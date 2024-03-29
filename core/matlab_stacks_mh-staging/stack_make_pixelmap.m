function tc = stack_make_pixelmap(stack, varargin)
%
% tc is a var of size [nCells,nTimepts]
%
%$Id: stack_make_pixelmap.m 79 2008-01-17 17:05:40Z histed $

uDefs = { 'FrameTimeMs', [], ...
          'DoHighpassFilter', true, ...
          'HighpassHz', 0.05, ...
          'DoHighpassRestoreMean', true, ...
          'FilterType', 'firls', ...
          'FirlsFilterOrder', 300, ...
          'BaselineFrames', [], ...
          'FilterInterpFrames', [], ...
          'DoDFOF', true, ...
          'FrameNsToAverage', [], ...
        };



uo = stropt2struct(stropt_defaults(uDefs, varargin));

%%

[nCells, nTimepts] = size(tcIn);

tc = double(tcIn);

% precompute some things
if uo.DoDFOF || uo.DoHighpassRestoreMean
    if isempty(uo.BaselineFrames)
        error('Must specify baseline frames');
    end
    baseVals = mean(tc(:,uo.BaselineFrames),2);
    baseCols = repmat(baseVals, [1 nTimepts]);
end

%% filtering
if uo.DoHighpassFilter
    if isempty(uo.FrameTimeMs)
        error('Must specifiy FrameTimeMs in order to filter');
    end
    if ~strcmp(uo.FilterType, 'interp') && ~isempty(uo.FilterInterpFrames)
        error('Cannot specify FilterInterpFrames for non-''interp'' filters');
    end


    switch uo.FilterType
      case 'butter'
        % butterworth
        Fs = 1000/uo.FrameTimeMs;
        Fstop = uo.HighpassHz / (Fs/2);
        N = 6;
        [b,a] = butter(N, [Fstop], 'high');
        tcF = filtfilt(b,a,tc')';
      case 'firls'
        % finite impulse reponse
        Fs = 1000/uo.FrameTimeMs;
        fpts = [0 uo.HighpassHz/2 uo.HighpassHz Fs/2] / (Fs/2);
        apts = [0 0 1 1];
        wpts = [1 1];
        N = uo.FirlsFilterOrder;
        [b,a] = firls(N, fpts, apts, wpts);
        tcF = filtfilt(b,a,tc')';
      case 'interp'
        % compute filtered baseline
        allLowF = interp_lowpass(tc, uo.FilterInterpFrames, ...
                                 'butter', 1000/uo.FrameTimeMs, ...
                                 uo.HighpassHz);
        % subtract this baseline
        tcF = tc - allLowF;
      otherwise 
        error('unknown filter type %s', uo.FilterType);
    end

    if uo.DoHighpassRestoreMean
        % add mean back in
        hpMean = mean(tcF(:,uo.BaselineFrames),2);
        tcF = tcF + baseCols - repmat(hpMean, [1 nTimepts]);
    end

    tc = tcF;
end

%% dF/F
if uo.DoDFOF
    tc = (tc - baseCols) ./ baseCols;
end

%% averaging
if ~isempty(uo.FrameNsToAverage)
    [nFramesInAverage nOutFr] = size(uo.FrameNsToAverage);
    avgTc = repmat(NaN, [nCells, nOutFr]);
    for iF = 1:nOutFr
        avgTc(:,iF) = mean(tc(:,uo.FrameNsToAverage(:,iF)),2);
    end
    tc = avgTc;
end


% old effort to interp filter: subset, butterworth, re-interpolate.
% didn't work because interpolation is hard
%
% $$$       case 'interp'
% $$$ 
% $$$         % select a subset of frames, lowpass, interpolate, subtract
% $$$         assert(~isempty(uo.FilterInterpFrames));
% $$$         assert(all(diff(uo.FilterInterpFrames) > 0), ...
% $$$                'FilterInterpFrames must be monotonically increasing');
% $$$         
% $$$ 
% $$$         tcPre = tc(uo.FilterInterpFrames);
% $$$         
% $$$         % lowpass, butter
% $$$         Fs = 1000/uo.FrameTimeMs;
% $$$         Fstop = uo.HighpassHz / (Fs/2);
% $$$         N = 6;
% $$$         [b,a] = butter(N, [Fstop], 'low');
% $$$         lowF = filtfilt(b,a,tcPre')';
% $$$ 
% $$$         % interp to all frames
% $$$         allLowF = interp1(uo.FilterInterpFrames, lowF, 1:nTimepts, 'pchip');
% $$$         keyboard
% $$$         % subtract
% $$$         tcF = tc - allLowF;
% $$$         
% $$$         %% test, interp first
% $$$         r = interp1(uo.FilterInterpFrames, tc(:,uo.FilterInterpFrames), ...
% $$$                     1:nTimepts, 'pchip');
% $$$         
% $$$         %% end
