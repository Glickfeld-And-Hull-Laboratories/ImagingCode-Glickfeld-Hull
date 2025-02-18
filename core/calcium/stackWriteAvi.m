function stackWriteAvi(stack, outName, fps, tCmap)
%STACKWRITEAVI Writes stack as AVI movie
% stackWriteAvi(stack, outName, fps, tCmap)
%
%$Id: stackWriteAvi.m 295 2008-07-14 19:27:23Z vincent $

if nargin < 3, fps = 7; end
if nargin < 4, tCmap = gray(256); end

[nRows,nCols,nFrames,nPlanes] = size(stack);

aviOpts = { ...
    'compression', 'none', ...
    'fps', fps }; 
switch nPlanes
  case 1
    aviOpts = { aviOpts{:}, ...
                'colormap', tCmap };
  case 3
    % no colormap
  otherwise
    error('Invalid size of image: nPlanes %d', nPlanes);
end

% open file
aviobj = avifile(outName, aviOpts{:});

% try
    % add frames to the file
    for iF = 1:nFrames
        aviobj = addframe(aviobj, squeeze(stack(:,:,iF,:)));
    end
% catch ME
%     aviobj = close(aviobj);
% end

aviobj = close(aviobj);

