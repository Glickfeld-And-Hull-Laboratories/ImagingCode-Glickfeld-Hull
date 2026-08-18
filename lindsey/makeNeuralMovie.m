function makeNeuralMovie(data, dirStim, outputFile, climLow, climHigh)
% makeNeuralMovie  Write an .mp4 of neural activity with a drifting Gabor inset.
%
% INPUTS:
%   data        - [height x width x nFrames] neural activity movie (~500x700)
%   dirStim     - [nTrials x 3] matrix:
%                   col 1: grating direction (degrees, 0-360)
%                   col 2: frame stimulus appears (1-indexed)
%                   col 3: frame stimulus disappears (1-indexed)
%   outputFile  - output filename string, e.g. 'neural_movie.mp4'
%   climLow     - lower bound for data normalization (e.g. prctile(data(:),1))
%   climHigh    - upper bound for data normalization (e.g. prctile(data(:),99))
%
% EXAMPLE USAGE:
%   makeNeuralMovie(data, dirStim, 'output.mp4', prctile(data(:),1), prctile(data(:),99))

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
frameRate       = 15;       % Hz  - acquisition / playback frame rate
insetSize       = 100;      % pixels - inset window is 100x100
gratingDiameter = 80;       % pixels - visible grating diameter inside window
gratingFreq     = 2;        % Hz  - temporal frequency of drifting grating
spatialCycle    = 10;       % pixels per spatial cycle
margin          = 10;       % pixels from the top-right corner of the frame

% Gabor envelope: sigma chosen so the Gaussian fits naturally inside the
% grating diameter (sigma = diameter / 6 gives ~3-sigma falloff at the edge)
gabor_sigma     = gratingDiameter / 6;

% -------------------------------------------------------------------------
% Setup
% -------------------------------------------------------------------------
[height, width, nFrames] = size(data);

% Normalize data to [0,1] using user-supplied scaling range, then clip
dataNorm = (double(data) - climLow) / (climHigh - climLow);
dataNorm = max(0, min(1, dataNorm));
dataUint8 = uint8(dataNorm * 255);

% Build a lookup: for each frame, which trial is active? (-1 = none)
trialAtFrame = -ones(1, nFrames);
for t = 1:size(dirStim, 1)
    onFrame  = dirStim(t, 2);
    offFrame = dirStim(t, 3);
    trialAtFrame(onFrame:offFrame) = t;
end

% Inset position in the upper-right corner (1-indexed)
colStart = width  - margin - insetSize + 1;
colEnd   = width  - margin;
rowStart = margin;
rowEnd   = margin + insetSize - 1;

% Pre-build spatial coordinate grids centred on the inset
[xx, yy] = meshgrid(1:insetSize, 1:insetSize);
cx = (insetSize + 1) / 2;
cy = (insetSize + 1) / 2;

% Distance from centre (used for both circular clip and Gabor envelope)
dist = sqrt((xx - cx).^2 + (yy - cy).^2);

% Hard circular mask: clips everything outside gratingDiameter/2
circleRadius = gratingDiameter / 2;
circleMask   = dist <= circleRadius;

% Gabor envelope: Gaussian windowed to zero outside the hard circle
% The envelope is computed only once; direction rotates the carrier, not the envelope
gabor_envelope           = exp(-(dist.^2) / (2 * gabor_sigma^2));
gabor_envelope(~circleMask) = 0;   % hard-clip outside grating circle

% Mean-luminance gray (mid-grey in [0,1])
grayLevel = 0.5;

% -------------------------------------------------------------------------
% Video writer
% -------------------------------------------------------------------------
v = VideoWriter(outputFile, 'MPEG-4');
v.FrameRate = frameRate;
open(v);

% -------------------------------------------------------------------------
% Frame loop
% -------------------------------------------------------------------------
for f = 1:nFrames

    % Neural activity frame -> uint8 RGB (grayscale replicated to 3 channels)
    gray  = dataUint8(:, :, f);
    frame = repmat(gray, [1, 1, 3]);

    trial = trialAtFrame(f);

    if trial == -1
        % Inter-trial interval: mean-luminance gray square
        inset = grayLevel * ones(insetSize, insetSize);
    else
        % Active trial: drifting Gabor stimulus
        direction_deg = dirStim(trial, 1);
        direction_rad = deg2rad(direction_deg);

        % Project (x,y) onto grating direction to get spatial phase
        spatialPhase = (2 * pi / spatialCycle) * ...
            ((xx - cx) * cos(direction_rad) + (yy - cy) * sin(direction_rad));

        % Temporal phase advances from frame 1 (continuous, no per-trial reset)
        temporalPhase = 2 * pi * gratingFreq * (f / frameRate);

        % Sinusoidal carrier: 100% contrast, range [-1, 1]
        carrier = sin(spatialPhase - temporalPhase);

        % Gabor = mean luminance + Gaussian-modulated carrier
        % gabor_envelope ranges [0,1]; carrier ranges [-1,1]
        % Result stays within [0,1] because envelope <= 1
        inset = grayLevel + grayLevel * (gabor_envelope .* carrier);

        % Fill area outside grating circle with mean luminance
        inset(~circleMask) = grayLevel;
    end

    % Convert inset to uint8 and composite into all 3 colour channels
    insetUint8 = uint8(inset * 255);
    frame(rowStart:rowEnd, colStart:colEnd, 1) = insetUint8;
    frame(rowStart:rowEnd, colStart:colEnd, 2) = insetUint8;
    frame(rowStart:rowEnd, colStart:colEnd, 3) = insetUint8;

    writeVideo(v, frame);
end

close(v);
fprintf('Saved %d frames to %s\n', nFrames, outputFile);
end
