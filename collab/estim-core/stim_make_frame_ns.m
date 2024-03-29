function os = stim_make_frame_ns(nTotalFrames, stimEvery, preStimNs, ...
                                 postStimNs, stimsToLump, skipRepsAtStart, ...
                                 keepOnlyNs)
                                 
%
%$Id: stim_make_frame_ns.m 320 2008-08-22 13:58:05Z histed $

if nargin < 5, stimsToLump = 1; end
if nargin < 6, skipRepsAtStart = 0; end
if nargin < 7, keepOnlyNs = []; end
    

os.nBaseFrames = stimEvery-1;
if length(stimEvery) == 1
    os.stimNs = stimEvery:stimEvery:(nTotalFrames-1);
else
    % can specify a vector
    assert(isvector(stimEvery));
    assert(stimsToLump == 1);
    os.stimNs = stimEvery;
end

assert(length(os.stimNs) > 1);
assert(max(os.stimNs(:)) < nTotalFrames, ...
       'some stimNs greater than total n frames');
if ~isempty(skipRepsAtStart) && skipRepsAtStart > 0
    os.stimNs = os.stimNs(skipRepsAtStart+1:end);
elseif ~isempty(keepOnlyNs)
    os.stimNs = os.stimNs(keepOnlyNs);
end


os.nStims = length(os.stimNs);

nOutLumps = os.nStims ./ stimsToLump;
assert(round(nOutLumps) == nOutLumps, ...
       'stimsToLump not a multiple of total nStims: are params correct?');
%fprintf(1, '%s: n reps to average is %d\n', ...
%        mfilename, nOutLumps);
os.blockStartNs = repmat(NaN, 1, nOutLumps);
os.blockEndNs = repmat(NaN, 1, nOutLumps);
os.blockFrameNs = {};
for iS = 1:nOutLumps
    tFirstStimN = os.stimNs( (iS-1)*(stimsToLump)+1);
    tLastStimN = os.stimNs(iS*(stimsToLump));
    os.blockStartNs(iS) = tFirstStimN - preStimNs;
    os.blockEndNs(iS) = tLastStimN + postStimNs;
    os.blockFrameNs{iS} = os.blockStartNs(iS):os.blockEndNs(iS);

    tFN = 1:length(os.blockFrameNs{iS});
    blockStimNs{iS} = tFN(ismember(os.blockFrameNs{iS}, os.stimNs));    
end

os.frameNsToAverage = cat(1,os.blockFrameNs{:});
[os.nFramesPerAvg,os.nFramesInBlock] = size(os.frameNsToAverage);

% process blockStimNs: make sure they are all the same and output
colDiffs = range(cat(1, blockStimNs{:}),1);
assert(all(colDiffs) == 0);
os.blockStimNs = blockStimNs{1};


os.preStimNs = preStimNs;
os.postStimNs = postStimNs;
os.nStimsInSer = stimsToLump;
% size [nFramesPerAvg,nFramesInBlock]


